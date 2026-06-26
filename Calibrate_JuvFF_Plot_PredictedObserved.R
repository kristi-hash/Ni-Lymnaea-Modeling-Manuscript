# PREDICTED VS OBSERVED LENGTH OVER TIME — FoodFactor Calibration
# Calibration and plotting use the Ni0 (control) treatment only.
# Sim Copper = 0 is matched to observed Treatment == "Ni0".
#
# Calibrates FoodFactor separately for each test:
#   FH   — single FoodFactor, from FreshlyHatched_CalFoodFactor
#   2WFF — two FoodFactors (FF1, FF2), from 2WeekOld_CalFoodFactor
#   2WL  — two FoodFactors (FF1, FF2), from 2WeekOld_CalFoodFactor
#
# 2W two-phase food factor — two-stage simulation design:
#   FF1 (days 1–14): calibrated from sim_file_2W (FF1 randomised) using a
#             unified calibration pooling both 2WFF and 2WL Ni0 day-14
#             observations (equal weighting). Fixed FF1 value carried forward.
#   FF2 (days 15+):  calibrated separately per test from sim_file_2W2, where
#             FF1 is already fixed at the calibrated value and only FF2 is
#             randomised. No conditioning step needed — every iteration starts
#             from the correct day-14 state.
#
# Two multipanel figures:
#   Fig 1: Per-test calibration — single best
#   Fig 2: Per-test calibration — median of best 100
# All predicted lines as dashed
# Unit conversions: Observed µm -> mm (/1000) | Sim already in mm
# Ribbon = observed min-max | Points = observed mean

rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)

#####################
# USER SETTINGS

sim_file_FH  <- "FreshlyHatched_CalFoodFactor"   # FH simulations only
sim_file_2W  <- "2WeekOld_CalFoodFactor"          # 2WFF and 2WL — FF1 randomised
sim_file_2W2 <- "2WeekOld_CalFoodFactor2"         # 2WFF and 2WL — FF1 fixed at calibrated value, FF2 randomised

obs_file <- "SnailLengths_Juveniles_SimAligned.csv"
n_best   <- 100

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
# Restrict to Ni0 (control) treatment only; sim has only Copper = 0.

obs <- read_csv(obs_file, show_col_types = FALSE)
obs <- obs[!is.na(obs$Size), ]
obs$Size_mm <- obs$Size / 1000

obs <- obs %>%
  filter(Treatment == "Ni0") %>%
  mutate(Panel = case_when(
    Test == "FH"   ~ "FreshlyHatched",
    Test == "2WFF" ~ "TwoWeek_FF",
    Test == "2WL"  ~ "TwoWeek_L",
    TRUE           ~ NA_character_
  )) %>%
  filter(!is.na(Panel))

cat("Observed data: Ni0 (control) treatment only\n")
cat("Panels and day ranges:\n")
for (g in c("FreshlyHatched", "TwoWeek_FF", "TwoWeek_L")) {
  days <- sort(unique(obs$Day[obs$Panel == g]))
  cat(" ", g, "— days:", paste(days, collapse = ", "), "\n")
}
cat("\n")

# obs_mean: used for SSE matching (one mean per Panel x Day)
obs_mean <- obs %>%
  group_by(Panel, Day) %>%
  summarise(obs_Size_mm = mean(Size_mm, na.rm = TRUE), .groups = "drop")

# obs_plot: mean + min/max ribbon, filtered to plotted day range
obs_plot <- obs %>%
  group_by(Panel, Day) %>%
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
# LOAD SIM FILES

load_sim <- function(file_path) {
  sim <- read_delim(file_path, delim = " ", show_col_types = FALSE)
  sim$PopulationLength <- as.numeric(sim$PopulationLength)
  sim <- sim[!is.na(sim$PopulationLength) & sim$PopulationLength > 0, ]
  sim
}

cat("Loading FH simulation file...\n")
sim_fh_raw <- load_sim(sim_file_FH)

cat("Loading 2W simulation file (FF1 randomised)...\n")
sim_2w_raw <- load_sim(sim_file_2W)

cat("Loading 2W simulation file (FF2 randomised, FF1 fixed)...\n")
sim_2w2_raw <- load_sim(sim_file_2W2)

#####################
# MATCH SIM TO OBSERVED (control only: Copper = 0 matched to Ni0 obs_mean)

match_sim_to_obs <- function(sim_raw, panel_name) {
  obs_panel <- obs_mean %>% filter(Panel == panel_name)
  
  sim_raw %>%
    filter(Copper == 0) %>%
    inner_join(obs_panel, by = "Day") %>%
    mutate(Panel = panel_name)
}

sim_FH   <- match_sim_to_obs(sim_fh_raw, "FreshlyHatched")
sim_2WFF <- match_sim_to_obs(sim_2w_raw, "TwoWeek_FF")
sim_2WL  <- match_sim_to_obs(sim_2w_raw, "TwoWeek_L")

for (panel in c("FreshlyHatched", "TwoWeek_FF", "TwoWeek_L")) {
  dat <- get(c(FreshlyHatched = "sim_FH", TwoWeek_FF = "sim_2WFF",
               TwoWeek_L = "sim_2WL")[[panel]])
  if (nrow(dat) == 0) stop(
    "No rows matched for panel ", panel,
    " — check sim Copper = 0 exists and observed Ni0 days align."
  )
  cat("Panel:", panel, "| matched rows:", nrow(dat),
      "| days:", paste(sort(unique(dat$Day)), collapse = ", "), "\n")
}
cat("\n")

#####################
# COMPUTE SSE
# FH:  single SSE over all matched days
# 2W:  SSE_FF1 — Day == 14 only  (calibrates FF1, pre-switch growth)
#      SSE_FF2 — Day >= 14       (calibrates FF2, exposure-period growth;
#                                  Day 14 included to anchor the switch point)

add_sse_fh <- function(sim_panel) {
  sim_panel %>%
    mutate(
      SSE   = (PopulationLength - obs_Size_mm)^2,
      SSEid = as.character(Iteration)
    )
}

add_sse_2w <- function(sim_panel) {
  sim_panel %>%
    mutate(
      SSE     = (PopulationLength - obs_Size_mm)^2,
      SSEid   = as.character(Iteration),
      SSE_FF1 = ifelse(Day == 14, SSE, NA_real_),
      SSE_FF2 = ifelse(Day >= 14, SSE, NA_real_)
    )
}

sim_FH   <- add_sse_fh(sim_FH)
sim_2WFF <- add_sse_2w(sim_2WFF)
sim_2WL  <- add_sse_2w(sim_2WL)

#####################
# CALIBRATION HELPERS

# FH: single FoodFactor — LL over all matched days
calibrate_fh <- function(sim_panel, label) {
  sse_vec  <- tapply(sim_panel$SSE, sim_panel$SSEid, sum)
  LL       <- (-1/2) * log(sse_vec)
  LL_order <- sort(LL, decreasing = TRUE)
  best_n   <- head(LL_order, n = n_best)
  top_iters <- as.integer(names(best_n))
  
  ff_lookup <- sim_fh_raw %>%
    distinct(Iteration, FoodFactor) %>%
    filter(Iteration %in% top_iters)
  ff_top <- ff_lookup$FoodFactor[match(top_iters, ff_lookup$Iteration)]
  
  best_iter   <- top_iters[1]
  median_iter <- top_iters[which.min(abs(ff_top - median(ff_top)))]
  
  list(
    label       = label,
    type        = "single",
    top_iters   = top_iters,
    iter_single = best_iter,
    iter_median = median_iter,
    FF_single   = sim_fh_raw$FoodFactor[sim_fh_raw$Iteration == best_iter][1],
    LL_single   = as.numeric(LL_order[1]),
    FF_median   = median(ff_top),
    FF_Q1       = quantile(ff_top, 0.025),
    FF_Q3       = quantile(ff_top, 0.975),
    LL_mean     = mean(best_n),
    LL_max      = max(best_n)
  )
}

# 2W: sequential conditional calibration with unified FF1
#
#   Step 1 — unified FF1: pool SSE_FF1 from both 2WFF and 2WL (equal weighting
#             via normalisation), calibrate on Day == 14. Top n_best iterations
#             are shared across both tests.
#
#   Step 2 — per-test FF2: for each test, calibrate FF2 on Day >= 14 restricted
#             to the shared top-FF1 iterations.

# Build ff_by_iter lookup (shared across both 2W tests)
ff_by_iter_2w <- sim_2w_raw %>%
  filter(Iteration > 0) %>%
  group_by(Iteration) %>%
  summarise(FF1 = FoodFactor[Day <= 14][1],
            FF2 = FoodFactor[Day >  14][1],
            .groups = "drop")

# Step 1: unified FF1 — pool normalised SSE_FF1 from both panels
norm_sse_col <- function(sim_panel, col) {
  mean_val <- mean(sim_panel[[col]], na.rm = TRUE)
  sim_panel[[paste0(col, "_norm")]] <- sim_panel[[col]] / mean_val
  sim_panel
}

sim_2WFF <- norm_sse_col(sim_2WFF, "SSE_FF1")
sim_2WL  <- norm_sse_col(sim_2WL,  "SSE_FF1")

sim_2w_pooled_ff1 <- bind_rows(sim_2WFF, sim_2WL)
sse_ff1_unified   <- tapply(sim_2w_pooled_ff1$SSE_FF1_norm,
                            sim_2w_pooled_ff1$SSEid, sum, na.rm = TRUE)
LL_ff1   <- (-1/2) * log(sse_ff1_unified)
LL_ff1_o <- sort(LL_ff1, decreasing = TRUE)
top_ff1  <- as.integer(names(head(LL_ff1_o, n_best)))

ff1_top_vals <- ff_by_iter_2w$FF1[match(top_ff1, ff_by_iter_2w$Iteration)]

cal_FF1_unified <- list(
  top_iters = top_ff1,
  FF_single = ff_by_iter_2w$FF1[match(top_ff1[1], ff_by_iter_2w$Iteration)],
  FF_median = median(ff1_top_vals, na.rm = TRUE),
  FF_Q1     = quantile(ff1_top_vals, 0.025, na.rm = TRUE),
  FF_Q3     = quantile(ff1_top_vals, 0.975, na.rm = TRUE),
  LL_single = as.numeric(LL_ff1_o[1]),
  LL_mean   = mean(head(LL_ff1_o, n_best)),
  LL_max    = max(head(LL_ff1_o, n_best))
)

cat("Unified 2W FF1 calibration — top", n_best, "iterations selected.\n")
cat("  FF1 single best:", round(cal_FF1_unified$FF_single, 6), "\n")
cat("  FF1 median:     ", round(cal_FF1_unified$FF_median, 6), "\n\n")

# Step 2: per-test FF2 — calibrated directly from sim_2w2_raw where FF1 is
# already fixed at the calibrated value. No conditioning on top-FF1 iterations
# needed; every iteration starts from the correct day-14 snail state.

# Match sim_2w2_raw to observed for each 2W panel (Day >= 14 only for FF2)
sim_2WFF2 <- match_sim_to_obs(sim_2w2_raw, "TwoWeek_FF")
sim_2WL2  <- match_sim_to_obs(sim_2w2_raw, "TwoWeek_L")
sim_2WFF2 <- add_sse_2w(sim_2WFF2)
sim_2WL2  <- add_sse_2w(sim_2WL2)

ff_by_iter_2w2 <- sim_2w2_raw %>%
  filter(Iteration > 0) %>%
  group_by(Iteration) %>%
  summarise(FF1 = FoodFactor[Day <= 14][1],
            FF2 = FoodFactor[Day >  14][1],
            .groups = "drop")

calibrate_ff2 <- function(sim_panel, label) {
  sse_ff2  <- tapply(sim_panel$SSE_FF2, sim_panel$SSEid, sum, na.rm = TRUE)
  LL_ff2   <- (-1/2) * log(sse_ff2)
  LL_ff2_o <- sort(LL_ff2, decreasing = TRUE)
  top_ff2  <- as.integer(names(head(LL_ff2_o, n_best)))
  
  ff2_top_vals <- ff_by_iter_2w2$FF2[match(top_ff2, ff_by_iter_2w2$Iteration)]
  
  best_iter   <- top_ff2[1]
  median_iter <- top_ff2[which.min(abs(ff2_top_vals - median(ff2_top_vals, na.rm = TRUE)))]
  
  list(
    label       = label,
    type        = "two_phase",
    best_iter   = best_iter,
    median_iter = median_iter,
    FF1         = cal_FF1_unified,
    FF2 = list(
      top_iters = top_ff2,
      FF_single = ff_by_iter_2w2$FF2[match(best_iter,   ff_by_iter_2w2$Iteration)],
      FF_median = median(ff2_top_vals, na.rm = TRUE),
      FF_Q1     = quantile(ff2_top_vals, 0.025, na.rm = TRUE),
      FF_Q3     = quantile(ff2_top_vals, 0.975, na.rm = TRUE),
      LL_single = as.numeric(LL_ff2_o[1]),
      LL_mean   = mean(head(LL_ff2_o, n_best)),
      LL_max    = max(head(LL_ff2_o, n_best))
    )
  )
}

cal_FH   <- calibrate_fh(sim_FH, "FreshlyHatched")
cal_2WFF <- calibrate_ff2(sim_2WFF2, "TwoWeek_FF")
cal_2WL  <- calibrate_ff2(sim_2WL2,  "TwoWeek_L")

cal_group <- list(
  FreshlyHatched = cal_FH,
  TwoWeek_FF     = cal_2WFF,
  TwoWeek_L      = cal_2WL
)

#####################
# HELPER: get predicted lines for a given panel and iteration

get_sim_lines_panel <- function(panel_name, iter, sim_raw) {
  sim_raw %>%
    filter(Iteration == iter, Copper == 0) %>%
    group_by(Day) %>%
    summarise(sim_mean = mean(PopulationLength, na.rm = TRUE), .groups = "drop") %>%
    mutate(Panel = panel_name) %>%
    rowwise() %>%
    filter(Day >= t_exposure_map[Panel],
           is.na(day_max_map[Panel]) | Day <= day_max_map[Panel]) %>%
    ungroup()
}

# Predicted lines — single best and median of top 100
sim_single <- bind_rows(
  get_sim_lines_panel("FreshlyHatched", cal_FH$iter_single,   sim_fh_raw),
  get_sim_lines_panel("TwoWeek_FF",     cal_2WFF$best_iter,   sim_2w2_raw),
  get_sim_lines_panel("TwoWeek_L",      cal_2WL$best_iter,    sim_2w2_raw)
)

sim_median <- bind_rows(
  get_sim_lines_panel("FreshlyHatched", cal_FH$iter_median,   sim_fh_raw),
  get_sim_lines_panel("TwoWeek_FF",     cal_2WFF$median_iter, sim_2w2_raw),
  get_sim_lines_panel("TwoWeek_L",      cal_2WL$median_iter,  sim_2w2_raw)
)

#####################
# PLOTTING FUNCTION
# No colour scale needed — single treatment (control) per panel.

make_panel <- function(panel_id, panel_title, obs_data, sim_data,
                       show_x_label = FALSE) {
  
  obs_p <- obs_data %>% filter(Panel == panel_id)
  sim_p <- sim_data %>% filter(Panel == panel_id)
  
  ggplot() +
    geom_ribbon(
      data  = obs_p,
      aes(x = Day, ymin = obs_min, ymax = obs_max),
      fill  = "steelblue", alpha = 0.20
    ) +
    geom_point(
      data  = obs_p,
      aes(x = Day, y = obs_mean),
      colour = "steelblue", size = 2.5, shape = 16
    ) +
    geom_line(
      data      = sim_p,
      aes(x = Day, y = sim_mean),
      colour    = "steelblue", linewidth = 0.7, linetype = "dashed"
    ) +
    labs(title = panel_title,
         x     = if (show_x_label) "Day" else NULL,
         y     = "Shell length (mm)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 11),
      axis.title.x = if (show_x_label) element_text() else element_blank()
    )
}

#####################
# ASSEMBLE FIGURES

make_figure <- function(sim_data, fig_title) {
  p1 <- make_panel("FreshlyHatched", "A  Freshly Hatched",
                   obs_plot, sim_data, show_x_label = FALSE)
  p2 <- make_panel("TwoWeek_FF", "B  Two Week — Fish Flakes",
                   obs_plot, sim_data, show_x_label = FALSE)
  p3 <- make_panel("TwoWeek_L", "C  Two Week — Lettuce",
                   obs_plot, sim_data, show_x_label = TRUE)
  
  fig <- plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr")
  ggdraw(fig) +
    draw_label(fig_title, x = 0.5,  y = 0.995, hjust = 0.5,
               vjust = 1, fontface = "bold", size = 12) +
    draw_label("Day",     x = 0.45, y = 0.01,  size = 11)
}

fig1 <- make_figure(sim_single,
                    "Figure 1 — Per-test calibration (best simulation)")
fig2 <- make_figure(sim_median,
                    "Figure 2 — Per-test calibration (median of best 100)")

print(fig1)
print(fig2)

#####################
# CONSOLE CALIBRATION REPORT

print_cal_fh <- function(cal) {
  cat("-------------------------------------------------------------\n")
  cat(" Test:", cal$label, "\n")
  cat("-------------------------------------------------------------\n")
  cat(" Single best (max LL):\n")
  cat("   FoodFactor:", round(cal$FF_single, 6), "\n")
  cat("   LL:        ", round(cal$LL_single,  4), "\n\n")
  cat(" Median of top", n_best, "fits:\n")
  cat("   FoodFactor:", round(cal$FF_median, 6),
      "  [2.5%:", round(cal$FF_Q1, 6), "— 97.5%:", round(cal$FF_Q3, 6), "]\n")
  cat("   Mean LL:   ", round(cal$LL_mean, 4), "\n")
  cat("   Max  LL:   ", round(cal$LL_max,  4), "\n\n")
}

print_cal_2w <- function(cal) {
  cat("-------------------------------------------------------------\n")
  cat(" Test:", cal$label, "\n")
  cat("-------------------------------------------------------------\n")
  cat(" FF2 calibration (Day >= 14, within top-FF1 iterations):\n")
  cat("   Single best — FoodFactor:", round(cal$FF2$FF_single, 6),
      "  LL:", round(cal$FF2$LL_single, 4), "\n")
  cat("   Median top", n_best, "—  FoodFactor:", round(cal$FF2$FF_median, 6),
      "  [2.5%:", round(cal$FF2$FF_Q1, 6), "— 97.5%:", round(cal$FF2$FF_Q3, 6), "]\n")
  cat("   Mean LL:", round(cal$FF2$LL_mean, 4),
      " | Max LL:", round(cal$FF2$LL_max, 4), "\n\n")
}

cat("\n=============================================================\n")
cat(" CALIBRATION REPORT — FoodFactor (Ni0 control only)\n")
cat("=============================================================\n\n")

cat("-------------------------------------------------------------\n")
cat(" FH — FoodFactor\n")
cat("-------------------------------------------------------------\n")
print_cal_fh(cal_FH)

cat("-------------------------------------------------------------\n")
cat(" 2W — Unified FF1 (pooled 2WFF + 2WL, Day == 14)\n")
cat("-------------------------------------------------------------\n")
cat(" Single best — FoodFactor:", round(cal_FF1_unified$FF_single, 6),
    "  LL:", round(cal_FF1_unified$LL_single, 4), "\n")
cat(" Median top", n_best, "—  FoodFactor:", round(cal_FF1_unified$FF_median, 6),
    "  [2.5%:", round(cal_FF1_unified$FF_Q1, 6),
    "— 97.5%:", round(cal_FF1_unified$FF_Q3, 6), "]\n")
cat(" Mean LL:", round(cal_FF1_unified$LL_mean, 4),
    " | Max LL:", round(cal_FF1_unified$LL_max, 4), "\n\n")

print_cal_2w(cal_2WFF)
print_cal_2w(cal_2WL)