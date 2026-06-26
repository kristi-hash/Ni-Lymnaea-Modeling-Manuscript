# PREDICTED VS OBSERVED LENGTH OVER TIME — Growth PMoA
# Four multipanel figures:
#   Fig 1: Combined calibration — single best
#   Fig 2: Combined calibration — median of best 100
#   Fig 3: Individual calibrations — single best
#   Fig 4: Individual calibrations — median of best 100
# All predicted lines as dashed
# Unit conversions: Observed µm -> mm (/1000) | Sim already in mm
# Ribbon = observed min-max | Points = observed mean
# Calibration option: recode lowest observed Ni to 0, matched to sim Ni=0
#   Lowest non-zero sim concentration excluded automatically (no matching obs)

rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)  # robust multipanel assembly (replaces plot_grid)
library(gridExtra)  # for PDF text page

#####################
# USER SETTINGS

sim_file  <- "Juveniles_CalibrateEP_FreeIon_Growth_A0to1_B0to5"
obs_file  <- "SnailLengths_Juveniles_SimAligned.csv"
n_best    <- 100

# Select which observed Ni column to calibrate against:
#   "DissNi"  — dissolved Ni concentration
#   "FreeIon" — free ion concentration
ni_source <- "FreeIon"

# Calibration matching option:
#   FALSE — standard: each observed Ni matched to same sim Copper value
#   TRUE  — recode lowest observed Ni per panel to 0, matched to sim Ni=0;
#            lowest non-zero sim concentration excluded (no matching obs);
#            observations and sim Ni=0 both plotted as 0 µg/L
match_lowest_to_control <- FALSE

# PDF output:
#   TRUE  — save all four figures and the calibration report to a PDF file
#   FALSE — print figures to the R graphics device only
save_pdf <- FALSE
pdf_file <- paste0("Calibration Result ", sim_file, ".pdf")

t_exposure_map <- c(
  FreshlyHatched = 0,
  TwoWeek_FF     = 14,
  TwoWeek_L      = 14
)

day_max_map <- c(
  FreshlyHatched = 28,
  TwoWeek_FF     = NA,
  TwoWeek_L      = NA
)

#####################
# LOAD AND PREPARE OBSERVED DATA

obs <- read_csv(obs_file, show_col_types = FALSE)
obs <- obs[!is.na(obs$Size), ]
obs$Size_mm <- obs$Size / 1000

if (!ni_source %in% colnames(obs)) stop(
  "Column '", ni_source, "' not found in ", obs_file,
  ".\nAvailable columns: ", paste(colnames(obs), collapse = ", ")
)

obs$Ni <- obs[[ni_source]]

obs <- obs %>%
  mutate(Panel = case_when(
    Test == "FH"   ~ "FreshlyHatched",
    Test == "2WFF" ~ "TwoWeek_FF",
    Test == "2WL"  ~ "TwoWeek_L",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Panel), Ni > 0, !is.na(Ni))

# Recode lowest Ni per panel to 0 for both calibration and plotting
if (match_lowest_to_control) {
  obs <- obs %>%
    group_by(Panel) %>%
    mutate(Ni = ifelse(Ni == min(Ni), 0, Ni)) %>%
    ungroup()
  cat("match_lowest_to_control = TRUE\n")
  cat("Lowest Ni per panel recoded to 0 for calibration and plotting.\n\n")
}

# Ni levels from data — includes 0 if recoded
ni_levels_map <- obs %>%
  group_by(Panel) %>%
  summarise(levels = list(sort(unique(Ni))), .groups = "drop") %>%
  { setNames(.$levels, .$Panel) }

cat("Using Ni source:", ni_source, "\n")
cat("Observed Ni levels for plotting:\n")
for (g in names(ni_levels_map)) {
  cat(" ", g, ":", ni_levels_map[[g]], "\n")
}
cat("\n")

obs_mean <- obs %>%
  group_by(Panel, Day, Ni) %>%
  summarise(obs_Size_mm = mean(Size_mm, na.rm = TRUE), .groups = "drop")

obs_plot <- obs %>%
  group_by(Panel, Day, Ni) %>%
  summarise(
    obs_mean = mean(Size_mm, na.rm = TRUE),
    obs_min  = min(Size_mm,  na.rm = TRUE),
    obs_max  = max(Size_mm,  na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  filter(Day >= t_exposure_map[Panel],
         is.na(day_max_map[Panel]) | Day <= day_max_map[Panel]) %>%
  ungroup()

#####################
# LOAD AND PREPARE SIMULATION OUTPUT
# BUG FIX: explicitly cast Copper, A_eq, B_eq to numeric after loading.
# read_delim reads these as character when BehaviorSpace repeats header rows;
# as.numeric() produces NA for those header-repeat rows, which are then dropped.

sim_raw_full <- read_delim(sim_file, delim = " ", show_col_types = FALSE) %>%
  mutate(
    PopulationLength = as.numeric(PopulationLength),
    Copper           = as.numeric(Copper),
    A_eq             = as.numeric(A_eq),
    B_eq             = as.numeric(B_eq),
    Iteration        = as.numeric(Iteration)
  ) %>%
  filter(!is.na(PopulationLength), PopulationLength > 0,
         !is.na(A_eq), !is.na(B_eq), !is.na(Copper)) %>%
  mutate(Panel = case_when(
    LifestageTest == "FH"                               ~ "FreshlyHatched",
    LifestageTest == "2WFF" | FoodType == "Fish flakes" ~ "TwoWeek_FF",
    LifestageTest == "2WL"  | FoodType == "Lettuce"     ~ "TwoWeek_L",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Panel))

# Match sim Copper to obs_mean Ni — unmatched sim concentrations drop naturally
sim <- sim_raw_full %>%
  left_join(obs_mean, by = c("Panel", "Day", "Copper" = "Ni")) %>%
  filter(!is.na(obs_Size_mm))

if (nrow(sim) == 0) stop(
  "No rows matched — check that sim Copper values align with ",
  ni_source, " values in the observed data."
)

cat("Matched rows:", nrow(sim), "\n")
cat("Sim Copper values matched:", sort(unique(sim$Copper)), "\n\n")

sim$SSE   <- (sim$PopulationLength - sim$obs_Size_mm)^2
sim$SSEid <- paste(sim$Iteration, sim$A_eq, sim$B_eq)

#####################
# CALIBRATION HELPER

calibrate <- function(sse_named_vec, label) {
  LL       <- (-1/2) * log(sse_named_vec)
  LL_order <- sort(LL, decreasing = TRUE)
  best_n   <- head(LL_order, n = n_best)
  
  sel <- strsplit(names(best_n), " ")
  df  <- data.frame(matrix(unlist(sel), nrow = length(best_n), byrow = TRUE),
                    stringsAsFactors = FALSE)
  colnames(df) <- c("Iteration", "A_eq", "B_eq")
  df$A_eq <- as.numeric(df$A_eq)
  df$B_eq <- as.numeric(df$B_eq)
  
  best_1      <- strsplit(names(LL_order)[1], " ")[[1]]
  A_eq_single <- as.numeric(best_1[2])
  B_eq_single <- as.numeric(best_1[3])
  
  list(
    label     = label,
    df_top100 = df,
    A_single  = A_eq_single,
    B_single  = B_eq_single,
    LL_single = as.numeric(LL_order[1]),
    A_median  = median(df$A_eq),
    A_Q1      = quantile(df$A_eq, 0.025),
    A_Q3      = quantile(df$A_eq, 0.975),
    B_median  = median(df$B_eq),
    B_Q1      = quantile(df$B_eq, 0.025),
    B_Q3      = quantile(df$B_eq, 0.975),
    LL_mean   = mean(best_n),
    LL_max    = max(best_n)
  )
}

#####################
# 1. OVERALL CALIBRATION (equal group weighting)

group_mean_SSE <- sim %>%
  group_by(Panel) %>%
  summarise(group_mean_SSE = mean(SSE, na.rm = TRUE), .groups = "drop")

sim <- sim %>%
  left_join(group_mean_SSE, by = "Panel") %>%
  mutate(SSE_norm = SSE / group_mean_SSE)

SSE_overall <- tapply(sim$SSE_norm, sim$SSEid, sum)
cal_overall <- calibrate(SSE_overall, "Overall (equal weighting)")

#####################
# 2. PER-GROUP CALIBRATION

groups    <- c("FreshlyHatched", "TwoWeek_FF", "TwoWeek_L")
cal_group <- lapply(groups, function(g) {
  sim_g <- sim[sim$Panel == g, ]
  sse_g <- tapply(sim_g$SSE, sim_g$SSEid, sum)
  calibrate(sse_g, g)
})
names(cal_group) <- groups

#####################
# HELPER: get predicted lines for a given A/B pair
# Uses sim_raw_full cached at load time — no re-read needed.

matched_sim_ni <- unique(obs_mean$Ni)

get_sim_lines <- function(A_val, B_val) {
  ab_grid  <- sim_raw_full %>%
    distinct(A_eq, B_eq) %>%
    mutate(dist = sqrt(((A_eq - A_val) / (abs(A_val) + 1e-12))^2 +
                         ((B_eq - B_val) / (abs(B_val) + 1e-12))^2))
  best_row <- ab_grid[which.min(ab_grid$dist), ]
  
  sim_raw_full %>%
    filter(A_eq == best_row$A_eq, B_eq == best_row$B_eq,
           Copper %in% matched_sim_ni) %>%
    group_by(Panel, Day, Copper) %>%
    summarise(sim_mean = mean(PopulationLength, na.rm = TRUE),
              .groups  = "drop") %>%
    rename(Ni = Copper) %>%
    rowwise() %>%
    filter(Day >= t_exposure_map[Panel],
           is.na(day_max_map[Panel]) | Day <= day_max_map[Panel]) %>%
    ungroup()
}

# Overall
sim_overall_single <- get_sim_lines(cal_overall$A_single, cal_overall$B_single)
sim_overall_median <- get_sim_lines(cal_overall$A_median, cal_overall$B_median)

# Per-group
sim_group_single <- bind_rows(lapply(groups, function(g) {
  get_sim_lines(cal_group[[g]]$A_single,
                cal_group[[g]]$B_single) %>% filter(Panel == g)
}))

sim_group_median <- bind_rows(lapply(groups, function(g) {
  get_sim_lines(cal_group[[g]]$A_median,
                cal_group[[g]]$B_median) %>% filter(Panel == g)
}))

#####################
# COLOUR SCALE
# Ni=0 labelled as "0 µg/L" when present

build_colour_scale <- function(panel_id) {
  ni_obs     <- ni_levels_map[[panel_id]]
  ni_labels  <- ifelse(ni_obs == 0, "0 µg/L", paste(ni_obs, "µg/L"))
  turbo_cols <- scales::viridis_pal(option = "turbo")(length(ni_obs))
  col_vals   <- setNames(turbo_cols, ni_labels)
  list(levels = ni_obs, labels = ni_labels, colours = col_vals)
}

#####################
# PLOTTING FUNCTION

make_panel <- function(panel_id, panel_title, obs_data, sim_data,
                       show_x_label = FALSE) {
  
  obs_p <- obs_data %>% filter(Panel == panel_id)
  sim_p <- sim_data %>% filter(Panel == panel_id)
  
  cs        <- build_colour_scale(panel_id)
  to_factor <- function(x) factor(x, levels = cs$levels, labels = cs$labels)
  
  obs_p$Ni_f <- to_factor(obs_p$Ni)
  sim_p$Ni_f <- to_factor(sim_p$Ni)
  
  obs_p <- obs_p %>% filter(!is.na(Ni_f))
  sim_p <- sim_p %>% filter(!is.na(Ni_f))
  
  legend_title <- if (ni_source == "DissNi") "Diss. Ni (µg/L)" else "Free Ion Ni (µg/L)"
  
  ggplot() +
    geom_ribbon(
      data  = obs_p,
      aes(x = Day, ymin = obs_min, ymax = obs_max, fill = Ni_f),
      alpha = 0.20
    ) +
    geom_point(
      data  = obs_p,
      aes(x = Day, y = obs_mean, colour = Ni_f),
      size  = 2.5, shape = 16
    ) +
    geom_line(
      data      = sim_p,
      aes(x = Day, y = sim_mean, colour = Ni_f),
      linewidth = 0.7, linetype = "dashed"
    ) +
    scale_colour_manual(name = legend_title, values = cs$colours,
                        drop = FALSE) +
    scale_fill_manual(  name = legend_title, values = cs$colours,
                        drop = FALSE) +
    labs(title = panel_title,
         x     = if (show_x_label) "Day" else NULL,
         y     = "Shell length (mm)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      axis.title.x    = if (show_x_label) element_text() else element_blank()
    )
}

#####################
# ASSEMBLE FIGURES
# BUG FIX (assembly): cowplot::plot_grid(align = "v", axis = "lr") performs
# gtable surgery to line up panel edges. When each sub-plot carries its own
# legend (as here — every panel has a different set of concentrations), the
# combined gtable can exceed the target device, collapsing the lower panels to
# ~zero height and clipping trailing legend rows. Symptom: only panel A drew,
# B and C blank, legend showed keys but lost most labels. Triggered on FreeIon
# (wider legend title + taller y-axis) but not DissNi on the same device.
#
# FIX: assemble with patchwork ( p1 / p2 / p3 ), which stacks and aligns axes
# without the collapse-prone surgery and keeps each panel's own legend.
# Earlier BUG FIX retained: no figure-level "Day" label — panel C carries the
# x-axis title via show_x_label = TRUE.

make_figure <- function(sim_data, fig_title) {
  p1 <- make_panel("FreshlyHatched", "A  Freshly Hatched",
                   obs_plot, sim_data, show_x_label = FALSE)
  p2 <- make_panel("TwoWeek_FF", "B  Two Week \u2014 Fish Flakes",
                   obs_plot, sim_data, show_x_label = FALSE)
  p3 <- make_panel("TwoWeek_L", "C  Two Week \u2014 Lettuce",
                   obs_plot, sim_data, show_x_label = TRUE)
  
  (p1 / p2 / p3) +
    plot_annotation(
      title = fig_title,
      theme = theme(plot.title = element_text(face = "bold", size = 12,
                                              hjust = 0.5))
    )
}

fig1 <- make_figure(sim_overall_single,
                    "Figure 1 \u2014 Combined calibration (best simulation)")
fig2 <- make_figure(sim_overall_median,
                    "Figure 2 \u2014 Combined calibration (median of best 100)")
fig3 <- make_figure(sim_group_single,
                    "Figure 3 \u2014 Individual calibrations (best simulation)")
fig4 <- make_figure(sim_group_median,
                    "Figure 4 \u2014 Individual calibrations (median of best 100)")

#####################
# CONSOLE CALIBRATION REPORT

print_cal <- function(cal) {
  cat("-------------------------------------------------------------\n")
  cat(" Group:", cal$label, "\n")
  cat("-------------------------------------------------------------\n")
  cat(" Single best (max LL):\n")
  cat("   A_eq:    ", round(cal$A_single, 6), "\n")
  cat("   B_eq:    ", round(cal$B_single, 4), "\n")
  cat("   LL:      ", round(cal$LL_single, 4), "\n\n")
  cat(" Median of top", n_best, "fits:\n")
  cat("   A_eq:    ", round(cal$A_median, 6),
      "  [2.5%:", round(cal$A_Q1, 6), "\u2014 97.5%:", round(cal$A_Q3, 6), "]\n")
  cat("   B_eq:    ", round(cal$B_median, 4),
      "  [2.5%:", round(cal$B_Q1, 4), "\u2014 97.5%:", round(cal$B_Q3, 4), "]\n")
  cat("   Mean LL: ", round(cal$LL_mean, 4), "\n")
  cat("   Max  LL: ", round(cal$LL_max,  4), "\n\n")
}

report_lines <- capture.output({
  cat("\n=============================================================\n")
  cat(" CALIBRATION REPORT \u2014 Growth PMoA\n")
  cat(" Ni source:", ni_source,
      "| match_lowest_to_control:", match_lowest_to_control, "\n")
  cat("=============================================================\n\n")
  print_cal(cal_overall)
  for (g in groups) print_cal(cal_group[[g]])
})

#####################
# OUTPUT — screen or PDF

if (save_pdf) {
  pdf(pdf_file, width = 8.5, height = 11)
  print(fig1)
  print(fig2)
  print(fig3)
  print(fig4)
  
  # Calibration report as a monospaced text page
  report_grob <- grid::textGrob(
    paste(report_lines, collapse = "\n"),
    x    = 0.05, y = 0.97,
    just = c("left", "top"),
    gp   = grid::gpar(fontfamily = "mono", fontsize = 7.5)
  )
  grid::grid.newpage()
  grid::grid.draw(report_grob)
  
  dev.off()
  cat("\nPDF saved to:", pdf_file, "\n")
} else {
  print(fig1)
  print(fig2)
  print(fig3)
  print(fig4)
}

cat(paste(report_lines, collapse = "\n"), "\n")