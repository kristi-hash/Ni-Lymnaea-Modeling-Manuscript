# CONTROL PREDICTIONS vs OBSERVED — Ni0, with individual variation
#
# Four panels: size thresholds > 5, > 15, > 25 mm (A–C) and cumulative reproduction (D),
# same layout/aesthetic. This is NOT a calibration: the simulation file is a single
# parameter set (MortalityConstant = 0.90, Ration = 9500, cv = 0.1) run as many
# stochastic replicates (iterations). Individual variation is shown as the spread
# across those replicates.
#
#   Observed (grey): mean = points joined by a DASHED line; ribbon = min–max across
#                    the Ni0 control tanks (6, 7, 14, 17, 19) at each timepoint.
#   Simulation (blue): mean = SOLID line; ribbon = min–max across replicates at each day.
#
# Panel D compares observed cumulative clutches (as-is) from Population_ReproductiveData.csv
# with the sim's Cumulativeoffspring / eggs_per_clutch (post-reset) as approximate clutches.
#
# Observed counts use Length >= X mm. The ~30 mm stocked snails clear the 25 mm
# threshold, and the sim reports Count25mm = 10 at Day 180, so observed and simulated
# panels both start at 10 (no boundary mismatch). 25 mm = OECD reproductive-adult minimum.
# Observed "Day" counts testing days; sim Day = observed Day + sim_day_offset (180).

rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#####################
# USER SETTINGS

sim_file <- "PopulationExperiment_ControlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts"
obs_file <- "PopulationDensity_PhotoLengths.csv"
rep_file <- "Population_ReproductiveData.csv"   # observed clutches (panel D)

# Panel D (cumulative reproduction), in CLUTCHES: observed cumulative clutches are
# used as-is; the sim's Cumulativeoffspring is divided by eggs_per_clutch to approximate
# clutches.
eggs_per_clutch <- 100
# The reproductive file dates are aligned to testing-day space via the stocking date
# (testing Day 0 = sim Day 180).
stock_date    <- as.Date("2023-05-15")
# Panel D anchor: the reproductive monitoring's own Day 0 is the first reading
# (2023-05-26), where cumulative clutches = 0. Anchor the OBSERVED reproduction series
# to this date so that zero reading sits at Day_test = 0 — matching the simulation,
# which begins accumulating from its own post-reset Day 0. (The file's "Day" column is
# not used because it drifts from true elapsed days later in the series; anchoring by
# date preserves correct spacing while putting the Day-0 reading at the origin.)
repro_start_date <- as.Date("2023-05-26")
# Cumulativeoffspring holds the 180-day burn-in total at sim Day 180 and is reset for
# the experiment phase, so the first valid post-reset value is sim Day 181. Exclude
# Day 180 (the burn-in total) — including it spikes the curve at the origin.
repro_reset_day <- 181

ni0_tanks    <- c(6, 7, 14, 17, 19)             # control replicate tanks
observedDays <- c(77, 91, 105, 119, 133)        # observed TESTING days (Day 0 has no lengths)

# Observed "Day" counts testing days only; the simulation's testing period begins
# at ExposureDay = 180 (sim Day 180 has CountIndividuals = 10, matching stocking).
# Work in TESTING-day space: Day_test = sim Day - sim_day_offset.
sim_day_offset <- 180

# Plot window (TESTING-day space), given in SIM-day terms for convenience. Observed
# points/ribbon only appear at observedDays; the simulation envelope spans the whole
# window. Starting at sim Day 180 (observation Day 0) includes the early reproduction
# surge, which will dominate the > 5 mm panel.
plot_sim_day_min <- 180                          # start at sim Day 180 (observation Day 0)
plot_sim_day_max <- 330                          # ... out to sim Day 330
plot_day_min <- plot_sim_day_min - sim_day_offset  # testing-day space (180 - 180 = 0)
plot_day_max <- plot_sim_day_max - sim_day_offset  # testing-day space (330 - 180 = 150)

panel_labels <- c(gt5 = "A  > 5 mm", gt15 = "B  > 15 mm", gt25 = "C  > 25 mm")
col_obs <- "grey40"
col_sim <- "blue"

save_pdf <- FALSE
pdf_file <- "ControlPredictions vs Observed Ni0.pdf"

# Day-0 stocking: each Ni0 tank was stocked with stock_n individuals at ~stock_length mm.
# Lengths weren't photographed at Day 0, so the known stocking counts are added directly.
stock_n      <- 10        # individuals stocked per tank at Day 0
stock_length <- 30        # approximate length (mm) of stocked individuals

#####################
# OBSERVED (Ni0 control): mean + min–max across tanks

obs <- read_csv(obs_file, show_col_types = FALSE)
names(obs) <- trimws(names(obs))

obs_plot <- obs %>%
  filter(Treatment == "Ni0", Tank %in% ni0_tanks, !is.na(Length)) %>%
  group_by(Day, Tank) %>%
  summarise(
    gt5  = sum(Length >= 5),
    gt15 = sum(Length >= 15),
    gt25 = sum(Length >= 25),
    .groups = "drop"
  ) %>%
  pivot_longer(c(gt5, gt15, gt25), names_to = "Thr", values_to = "count") %>%
  group_by(Thr, Day) %>%
  summarise(
    obs_mean = mean(count),
    obs_min  = min(count),
    obs_max  = max(count),
    .groups  = "drop"
  )

# Inject the Day-0 stocking observation: every stocked snail clears any threshold
# at or below stock_length, so each panel gets stock_n at Day 0 (identical across tanks).
thr_vals <- c(gt5 = 5, gt15 = 15, gt25 = 25)
day0_obs <- data.frame(
  Thr      = names(thr_vals),
  Day      = 0,
  obs_mean = ifelse(stock_length >= thr_vals, stock_n, 0),
  obs_min  = ifelse(stock_length >= thr_vals, stock_n, 0),
  obs_max  = ifelse(stock_length >= thr_vals, stock_n, 0),
  stringsAsFactors = FALSE,
  row.names = NULL
)
obs_plot <- bind_rows(obs_plot, day0_obs) %>% arrange(Thr, Day)

#####################
# OBSERVED REPRODUCTION (Ni0): cumulative clutches, mean + min–max across tanks
#
# NOTE: read from CSV (was .xlsx). The Date column is M/D/Y (e.g. 5/26/2023), so
# it is parsed with format = "%m/%d/%Y". The clutches column is `Cummulative Clutches`
# (double-m, with a space), matching the CSV header. Earliest date 5/26/2023 gives
# Day_test = 11 against stock_date 2023-05-15, consistent with the file's Day column
# starting 11 d after stocking.

rep_raw <- read_csv(rep_file, show_col_types = FALSE)
names(rep_raw) <- trimws(names(rep_raw))

rep_plot <- rep_raw %>%
  filter(Treatment == "Ni0", Tank %in% ni0_tanks) %>%
  mutate(
    Day_test = as.numeric(as.Date(Date, format = "%m/%d/%Y") - repro_start_date),  # repro Day 0 = first reading (5/26)
    clutches = `Cummulative Clutches`                          # observed clutches, as-is
  ) %>%
  group_by(Day_test) %>%
  summarise(
    obs_mean = mean(clutches),
    obs_min  = min(clutches),
    obs_max  = max(clutches),
    .groups  = "drop"
  )

#####################
# SIMULATION ENSEMBLE: mean + min–max across replicates (individual variation)
# Space-delimited; if the leading space mis-parses columns, use read_table(sim_file).

sim_raw <- read_delim(sim_file, delim = " ", show_col_types = FALSE) %>%
  mutate(
    Iteration = as.numeric(Iteration),
    Day       = as.numeric(Day),
    Count5mm  = as.numeric(Count5mm),
    Count15mm = as.numeric(Count15mm),
    Count25mm = as.numeric(Count25mm)
  ) %>%
  filter(!is.na(Iteration), Day > 0) %>%          # drop header-repeats and Day-0 seed block
  mutate(Day_test = Day - sim_day_offset)

sim_plot <- sim_raw %>%
  select(Iteration, Day_test, Count5mm, Count15mm, Count25mm) %>%
  pivot_longer(c(Count5mm, Count15mm, Count25mm),
               names_to = "CountCol", values_to = "sim") %>%
  mutate(Thr = recode(CountCol,
                      Count5mm = "gt5", Count15mm = "gt15", Count25mm = "gt25")) %>%
  filter(Day_test >= plot_day_min, Day_test <= plot_day_max) %>%
  group_by(Thr, Day_test) %>%
  summarise(
    sim_mean = mean(sim, na.rm = TRUE),
    sim_min  = min(sim,  na.rm = TRUE),
    sim_max  = max(sim,  na.rm = TRUE),
    n_rep    = dplyr::n(),
    .groups  = "drop"
  )

# Simulated cumulative reproduction (panel D): Cumulativeoffspring -> clutches
# (divide by eggs_per_clutch), post-reset only.
sim_repro <- sim_raw %>%
  mutate(sim_clutches = as.numeric(Cumulativeoffspring) / eggs_per_clutch) %>%
  filter(Day >= repro_reset_day, Day_test <= plot_day_max) %>%
  group_by(Day_test) %>%
  summarise(
    sim_mean = mean(sim_clutches, na.rm = TRUE),
    sim_min  = min(sim_clutches,  na.rm = TRUE),
    sim_max  = max(sim_clutches,  na.rm = TRUE),
    .groups  = "drop"
  )

#####################
# PLOTTING

make_panel <- function(thr_id, show_x = FALSE) {
  o <- obs_plot %>% filter(Thr == thr_id)
  s <- sim_plot %>% filter(Thr == thr_id)
  
  ggplot() +
    # Simulation: blue ribbon + solid mean line
    geom_ribbon(data = s, aes(x = Day_test, ymin = sim_min, ymax = sim_max),
                fill = col_sim, alpha = 0.15) +
    geom_line(data = s, aes(x = Day_test, y = sim_mean,
                            colour = "Simulation", linetype = "Simulation"),
              linewidth = 0.8) +
    # Observed: grey ribbon + dashed mean line + points
    geom_ribbon(data = o, aes(x = Day, ymin = obs_min, ymax = obs_max),
                fill = "grey70", alpha = 0.35) +
    geom_line(data = o, aes(x = Day, y = obs_mean,
                            colour = "Observed (Ni0)", linetype = "Observed (Ni0)"),
              linewidth = 0.7) +
    geom_point(data = o, aes(x = Day, y = obs_mean, colour = "Observed (Ni0)"),
               size = 2.2) +
    scale_colour_manual(name = NULL,
                        values = c("Observed (Ni0)" = col_obs, "Simulation" = col_sim)) +
    scale_linetype_manual(name = NULL,
                          values = c("Observed (Ni0)" = "dashed", "Simulation" = "solid")) +
    labs(title = panel_labels[thr_id],
         x = if (show_x) "Day" else NULL,
         y = "Number of snails") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "right",
      axis.title.x    = if (show_x) element_text() else element_blank()
    )
}

# Reproduction panel: same grey-observed / blue-simulation aesthetic, clutches on the y-axis.
make_repro_panel <- function(show_x = TRUE) {
  ggplot() +
    geom_ribbon(data = sim_repro, aes(x = Day_test, ymin = sim_min, ymax = sim_max),
                fill = col_sim, alpha = 0.15) +
    geom_line(data = sim_repro, aes(x = Day_test, y = sim_mean,
                                    colour = "Simulation", linetype = "Simulation"),
              linewidth = 0.8) +
    geom_ribbon(data = rep_plot, aes(x = Day_test, ymin = obs_min, ymax = obs_max),
                fill = "grey70", alpha = 0.35) +
    geom_line(data = rep_plot, aes(x = Day_test, y = obs_mean,
                                   colour = "Observed (Ni0)", linetype = "Observed (Ni0)"),
              linewidth = 0.7) +
    geom_point(data = rep_plot, aes(x = Day_test, y = obs_mean, colour = "Observed (Ni0)"),
               size = 2.2) +
    scale_colour_manual(name = NULL,
                        values = c("Observed (Ni0)" = col_obs, "Simulation" = col_sim)) +
    scale_linetype_manual(name = NULL,
                          values = c("Observed (Ni0)" = "dashed", "Simulation" = "solid")) +
    labs(title = paste0("D  Cumulative reproduction (sim offspring / ", eggs_per_clutch, ")"),
         x = if (show_x) "Day" else NULL,
         y = "Cumulative clutches") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "right",
      axis.title.x    = if (show_x) element_text() else element_blank()
    )
}

fig <- (make_panel("gt5",  show_x = FALSE) /
          make_panel("gt15", show_x = FALSE) /
          make_panel("gt25", show_x = FALSE) /
          make_repro_panel(show_x = TRUE)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Control predictions vs observed (Ni0) \u2014 individual variation",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  )

#####################
# OUTPUT — screen or PDF

if (save_pdf) {
  pdf(pdf_file, width = 8.5, height = 13)
  print(fig)
  dev.off()
  cat("PDF saved to:", pdf_file, "\n")
} else {
  print(fig)
}

# Quick check: replicate count contributing to each day (should be constant within
# the window if all replicates reach it).
cat("\nReplicates per day (window):\n")
print(sim_plot %>% group_by(Day_test) %>% summarise(n_rep = max(n_rep), .groups = "drop"))