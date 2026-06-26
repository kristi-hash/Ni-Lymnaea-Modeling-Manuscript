# EXPOSURE PREDICTIONS vs OBSERVED — cumulative reproduction, 6-panel grid
#
# Produces TWO 6-panel figures in one run: one for FreeIon, one for DissNi.
# Each figure is a 3 x 2 grid:
#   rows = treatments (Ni1, Ni2, Ni3) | LEFT column = OBSERVED | RIGHT column = SIMULATION
#
# Metric: CUMULATIVE CLUTCHES over time.
#   - Observed  : `Cummulative Clutches` per tank (Population_ReproductiveData.csv).
#   - Simulated : Cumulativeoffspring / eggs_per_clutch (approximate clutches).
# In every panel: control (grey) + relevant concentration (colour), each a mean SOLID
# line + min–max ribbon. Observed panels also show individual treatment tanks as points.
#
# Free ion and dissolved Ni are SEPARATE simulation runs; ni_sources selects which sim
# file(s) to read and which per-treatment concentrations to match against `Copper`.
# Observations are matched by treatment, so the observed column is identical across both.
#
# TIME ALIGNMENT (differs from the size scripts):
#   - Observed reproduction is anchored to its first reading (repro_start_date, 2023-05-26),
#     where cumulative clutches = 0; Day = days since that date.
#   - Simulated Cumulativeoffspring holds the 180-day burn-in total at sim Day 180 and
#     resets for the experiment phase, so the first valid post-reset day is repro_reset_day
#     (181). Day_test = sim Day - sim_day_offset (180), so the sim curve starts at Day 1.
# Incomplete (still-running) sim replicates are dropped; observed ribbons span the five
# replicate tanks. NOTE: simulated clutches (whole-population offspring/eggs_per_clutch)
# are typically on a larger scale than observed per-tank clutches — panels use independent
# y-axes so each is readable; compare shape/dose-response rather than absolute height.

rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

#####################
# USER SETTINGS

ni_sources <- c("FreeIon", "DissNi")    # produce a 6-panel figure for each

sim_files <- c(
  FreeIon = "PopulationEffectlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts_FreeIon_Growth",
  DissNi  = "PopulationEffectlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts_DissNi_Growth"
)
obs_file <- "Population_ReproductiveData.csv"

sim_day_offset  <- 180                 # sim Day 180 = testing Day 0 (stocking)
eggs_per_clutch <- 100                 # sim Cumulativeoffspring / this -> approx clutches
repro_reset_day <- 181                 # first valid post-reset sim Day (180 = burn-in total)
repro_start_date <- as.Date("2023-05-26")  # observed repro Day 0 (first reading, cum = 0)
control_conc    <- 0                   # control simulation Copper value (verify for DissNi)

# Treatment <-> free-ion conc <-> dissolved-Ni conc <-> colour.
exposures <- tibble(
  Treatment    = c("Ni1",      "Ni2",      "Ni3"),
  free_ion     = c(0.409,      1.342,      12.114),    # matches Copper in the FreeIon sim
  dissolved_ni = c(2.30,       7.13,       65.24),     # matches Copper in the DissNi sim
  colour       = c("#1f78b4",  "#ef8a06",  "#d7191c")
)

conc_label_metric <- "source"          # "source" | "free_ion" | "dissolved" | "both"
conc_units        <- "\u00b5g/L"

control_label  <- "Ni0 (control)"
control_colour <- "grey40"
control_fill   <- "grey70"

plot_day_min <- 0       # testing-day window (sim Day 180 -> 330)
plot_day_max <- 150

save_pdf <- FALSE
pdf_file <- "ExposurePredictions vs Observed Repro 6panel.pdf"

# Shared colour/fill scales so all panels feed ONE collected legend.
series_levels <- c(control_label, exposures$Treatment)
series_cols   <- setNames(c(control_colour, exposures$colour), series_levels)
series_fills  <- setNames(c(control_fill,   exposures$colour), series_levels)

#####################
# OBSERVED (computed once — matched by treatment; cumulative clutches per tank)

obs <- read_csv(obs_file, show_col_types = FALSE)
names(obs) <- trimws(names(obs))

all_treats <- c("Ni0", exposures$Treatment)

obs_tank <- obs %>%
  filter(Treatment %in% all_treats, !is.na(`Cummulative Clutches`)) %>%
  mutate(
    Day      = as.numeric(as.Date(Date, format = "%m/%d/%Y") - repro_start_date),
    clutches = `Cummulative Clutches`
  ) %>%
  select(Treatment, Day, Tank, clutches)

obs_plot <- obs_tank %>%
  group_by(Treatment, Day) %>%
  summarise(
    obs_mean = mean(clutches),
    obs_min  = min(clutches),
    obs_max  = max(clutches),
    .groups  = "drop"
  )

#####################
# HELPERS

conc_for <- function(src) if (src == "FreeIon") exposures$free_ion else exposures$dissolved_ni

conc_text <- function(i, src) {
  fi <- exposures$free_ion[i]; dn <- exposures$dissolved_ni[i]
  metric <- if (conc_label_metric == "source")
    (if (src == "FreeIon") "free_ion" else "dissolved")
  else conc_label_metric
  switch(metric,
         free_ion  = paste0("free ion = ", fi, " ", conc_units),
         dissolved = paste0("dissolved Ni = ", dn, " ", conc_units),
         both      = paste0("dissolved Ni = ", dn, " ", conc_units,
                            ", free ion = ", fi, " ", conc_units),
         stop("conc_label_metric must be 'source', 'dissolved', 'free_ion', or 'both'")
  )
}

base_theme <- function(show_x, show_y) {
  theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "right",
      axis.title.x    = if (show_x) element_text() else element_blank(),
      axis.title.y    = if (show_y) element_text() else element_blank()
    )
}

common_scales <- list(
  scale_colour_manual(name = NULL, values = series_cols,  limits = series_levels, drop = FALSE),
  scale_fill_manual(values = series_fills, limits = series_levels, drop = FALSE),
  guides(fill = "none"),
  coord_cartesian(xlim = c(plot_day_min, plot_day_max))
)

# LEFT column: observed control (Ni0) + observed treatment (+ tank points).
make_obs_panel <- function(trt, conc_lab, show_x = FALSE, show_y = TRUE) {
  ctrl <- obs_plot %>% filter(Treatment == "Ni0")
  o    <- obs_plot %>% filter(Treatment == trt)
  pts  <- obs_tank %>% filter(Treatment == trt)
  ggplot() +
    geom_ribbon(data = ctrl, aes(Day, ymin = obs_min, ymax = obs_max, fill = control_label),
                alpha = 0.22) +
    geom_ribbon(data = o,    aes(Day, ymin = obs_min, ymax = obs_max, fill = trt),
                alpha = 0.22) +
    geom_line(data = ctrl, aes(Day, obs_mean, colour = control_label), linewidth = 0.8) +
    geom_line(data = o,    aes(Day, obs_mean, colour = trt),           linewidth = 0.8) +
    geom_point(data = pts, aes(Day, clutches, colour = trt), size = 1.5, alpha = 0.7,
               position = position_jitter(width = 1.2, height = 0, seed = 1)) +
    common_scales +
    labs(title = paste0(trt, " \u2014 observed (", conc_lab, ")"),
         x = if (show_x) "Day" else NULL,
         y = "Cumulative clutches") +
    base_theme(show_x, show_y)
}

# RIGHT column: simulated control (Copper 0) + simulated matched exposure.
make_sim_panel <- function(trt, conc, conc_lab, sim_sum, show_x = FALSE, show_y = FALSE) {
  ctrl <- sim_sum %>% filter(Copper == control_conc)
  expo <- sim_sum %>% filter(Copper == conc)
  ggplot() +
    geom_ribbon(data = ctrl, aes(Day_test, ymin = sim_min, ymax = sim_max, fill = control_label),
                alpha = 0.22) +
    geom_ribbon(data = expo, aes(Day_test, ymin = sim_min, ymax = sim_max, fill = trt),
                alpha = 0.22) +
    geom_line(data = ctrl, aes(Day_test, sim_mean, colour = control_label), linewidth = 0.8) +
    geom_line(data = expo, aes(Day_test, sim_mean, colour = trt),           linewidth = 0.8) +
    common_scales +
    labs(title = paste0(trt, " \u2014 simulation (", conc_lab, ")"),
         x = if (show_x) "Day" else NULL,
         y = NULL) +
    base_theme(show_x, show_y)
}

# Per-concentration mean + min–max simulated cumulative clutches for one basis.
build_sim_sum <- function(src) {
  sim_file <- sim_files[[src]]
  conc     <- conc_for(src)
  
  sim_raw <- read_delim(sim_file, delim = " ", show_col_types = FALSE) %>%
    mutate(
      Iteration          = as.numeric(Iteration),
      Copper             = as.numeric(Copper),
      Day                = as.numeric(Day),
      Cumulativeoffspring = as.numeric(Cumulativeoffspring)
    ) %>%
    filter(!is.na(Iteration), Day > 0) %>%
    mutate(Day_test = Day - sim_day_offset)
  
  avail_conc <- sort(unique(sim_raw$Copper))
  cat("[", src, "] sim file:", sim_file, "\n")
  cat("[", src, "] concentrations available (Copper):",
      paste(avail_conc, collapse = ", "), "\n")
  sel_conc <- c(control_conc, conc)
  missing_conc <- sel_conc[!sel_conc %in% avail_conc]
  if (length(missing_conc) > 0)
    warning("[", src, "] selected concentration(s) not found: ",
            paste(missing_conc, collapse = ", "),
            " | available: ", paste(avail_conc, collapse = ", "))
  
  # Drop still-running (incomplete) replicates (group by Copper AND Iteration).
  sim_full_day <- max(sim_raw$Day, na.rm = TRUE)
  sim_raw <- sim_raw %>%
    group_by(Copper, Iteration) %>%
    filter(max(Day) >= sim_full_day) %>%
    ungroup()
  
  used_conc <- c(control_conc, conc)
  sim_sum <- sim_raw %>%
    mutate(clutches = Cumulativeoffspring / eggs_per_clutch) %>%
    filter(Copper %in% used_conc, Day >= repro_reset_day, Day_test <= plot_day_max) %>%
    group_by(Copper, Day_test) %>%
    summarise(
      sim_mean = mean(clutches, na.rm = TRUE),
      sim_min  = min(clutches,  na.rm = TRUE),
      sim_max  = max(clutches,  na.rm = TRUE),
      n_rep    = dplyr::n(),
      .groups  = "drop"
    )
  
  cat("[", src, "] complete replicates per concentration:\n")
  print(sim_sum %>% group_by(Copper) %>% summarise(reps = max(n_rep), .groups = "drop"))
  sim_sum
}

# Assemble the full 6-panel figure for one basis.
build_figure <- function(src) {
  sim_sum <- build_sim_sum(src)
  
  n_rows <- nrow(exposures)
  conc   <- conc_for(src)
  panels <- list()
  for (i in seq_len(n_rows)) {
    is_bottom <- (i == n_rows)
    clab <- conc_text(i, src)
    panels[[length(panels) + 1]] <- make_obs_panel(exposures$Treatment[i], clab,
                                                   show_x = is_bottom, show_y = TRUE)
    panels[[length(panels) + 1]] <- make_sim_panel(exposures$Treatment[i], conc[i],
                                                   clab, sim_sum,
                                                   show_x = is_bottom, show_y = FALSE)
  }
  
  wrap_plots(panels, ncol = 2, byrow = TRUE) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Exposure predictions vs observed \u2014 cumulative reproduction",
      subtitle = paste0("Concentration basis: ",
                        if (src == "FreeIon") "free ion" else "dissolved Ni"),
      theme = theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
                    plot.subtitle = element_text(size = 10, hjust = 0.5))
    )
}

#####################
# BUILD BOTH FIGURES (skip a basis whose sim file is missing/unreadable)

bad_src <- setdiff(ni_sources, names(sim_files))
if (length(bad_src) > 0) stop("ni_sources not in sim_files: ", paste(bad_src, collapse = ", "))

figs <- list()
for (src in ni_sources) {
  f <- tryCatch(build_figure(src),
                error = function(e) { warning("Skipping ", src, ": ", conditionMessage(e)); NULL })
  if (!is.null(f)) figs[[src]] <- f
}
if (length(figs) == 0) stop("No figures produced — check the sim files exist.")

#####################
# OUTPUT

if (save_pdf) {
  pdf(pdf_file, width = 10, height = 11)
  for (f in figs) print(f)
  dev.off()
  cat("\nPDF saved to:", pdf_file, "  (", length(figs), "figure(s) )\n")
} else {
  for (f in figs) print(f)
}