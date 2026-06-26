# TREATMENT EFFECT AS % OF CONTROL OVER TIME — individuals > 5 mm
#
# Produces TWO 3-panel figures in one run: one for FreeIon, one for DissNi.
# Each figure is split into one panel per treatment (Ni1, Ni2, Ni3). For each treatment
# and timepoint the metric (snails > 5 mm) is expressed as a PERCENTAGE OF THE MATCHED
# CONTROL:
#       value(treatment) / value(control) * 100
#   - Observed  (SOLID line + circles)   : mean across treatment tanks / mean across Ni0 tanks
#   - Predicted (DASHED line + triangles): mean across exposure reps   / mean across Cu=0 reps
# 100% = same as control (dotted reference line); < 100% reduction; > 100% increase.
#
# Free ion and dissolved Ni are SEPARATE simulation runs (the Copper column is reported
# in the respective metric), so each figure reads its own sim file and matches its own
# per-treatment concentrations (free_ion vs dissolved_ni). Observations are matched by
# treatment, so the observed side is identical across both figures.
#
# Treatment colours match the 6-panel figure (Ni1 blue, Ni2 orange, Ni3 red). Both series
# are sampled ONLY at the observation days (points), with straight lines connecting them;
# the x-axis ends at the last observation day. Day-0 stocking anchors every treatment at
# 100%. Incomplete (still-running) sim replicates are dropped.

rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(tibble)

#####################
# USER SETTINGS

# Which bases to render (one 3-panel figure each). Keep both to produce both at once.
ni_sources <- c("FreeIon", "DissNi")

sim_files <- c(
  FreeIon = "PopulationEffectlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts_FreeIon_Growth",
  DissNi  = "PopulationEffectlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts_DissNi_Growth"
)
obs_file <- "PopulationDensity_PhotoLengths.csv"

threshold_mm   <- 5     # size threshold (uses the sim's Count5mm column)
sim_day_offset <- 180   # sim Day 180 = testing Day 0 (stocking)
control_conc   <- 0     # Copper value of the control simulation (verify for DissNi file)

# Treatment <-> free-ion conc <-> dissolved-Ni conc <-> colour.
# The Copper values to match depend on the basis: free_ion for the FreeIon sim,
# dissolved_ni for the DissNi sim.
exposures <- tibble(
  Treatment    = c("Ni1",      "Ni2",      "Ni3"),
  free_ion     = c(0.409,      1.342,      12.114),   # matches Copper in the FreeIon sim
  dissolved_ni = c(2.30,       7.13,       65.24),    # matches Copper in the DissNi sim
  colour       = c("#1f78b4",  "#ef8a06",  "#d7191c")
)
control_treatment <- "Ni0"

# Concentration shown in panel (facet) labels: "source" (match the figure's basis),
# "free_ion", "dissolved", or "both".
conc_label_metric <- "source"
conc_units        <- "\u00b5g/L"

stock_n      <- 10      # individuals stocked per tank at Day 0
stock_length <- 30      # approximate stocked length (mm)
include_day0_anchor <- TRUE   # anchor observed lines at 100% on Day 0 (stocking)

plot_day_min <- 0       # testing-day window (sim Day 180 -> 330)
plot_day_max <- 150

# Per-treatment panels: FALSE = shared y-axis (effect magnitudes comparable across
# treatments); TRUE = independent y per panel (each trajectory maximally readable).
free_y_panels <- FALSE

save_pdf <- FALSE
pdf_file <- "TreatmentEffect PctOfControl gt5mm.pdf"

#####################
# OBSERVED (computed once — identical across bases; matched by treatment)

obs <- read_csv(obs_file, show_col_types = FALSE)
names(obs) <- trimws(names(obs))

all_treats <- c(control_treatment, exposures$Treatment)

obs_day <- obs %>%
  filter(!is.na(Length), Treatment %in% all_treats) %>%
  group_by(Treatment, Day, Tank) %>%
  summarise(n_gt = sum(Length >= threshold_mm), .groups = "drop") %>%
  group_by(Treatment, Day) %>%
  summarise(mean_count = mean(n_gt), .groups = "drop")

if (include_day0_anchor) {
  obs_day <- bind_rows(obs_day, tibble(
    Treatment  = all_treats,
    Day        = 0,
    mean_count = ifelse(stock_length >= threshold_mm, stock_n, 0)
  ))
}

# Observation days — both series are sampled at exactly these days.
obs_meas_days <- sort(unique(obs_day$Day))

ctrl_obs <- obs_day %>%
  filter(Treatment == control_treatment) %>%
  select(Day, ctrl = mean_count)

obs_pct <- obs_day %>%
  filter(Treatment %in% exposures$Treatment) %>%
  left_join(ctrl_obs, by = "Day") %>%
  mutate(pct = ifelse(!is.na(ctrl) & ctrl > 0, 100 * mean_count / ctrl, NA_real_),
         Source = "Observed") %>%
  filter(!is.na(pct)) %>%
  select(Treatment, Day, pct, Source)

#####################
# HELPERS

trt_cols <- setNames(exposures$colour, exposures$Treatment)

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

# PREDICTED % of control for one basis.
build_sim_pct <- function(src) {
  sim_file <- sim_files[[src]]
  conc     <- conc_for(src)
  
  sim_raw <- read_delim(sim_file, delim = " ", show_col_types = FALSE) %>%
    mutate(
      Iteration = as.numeric(Iteration),
      Copper    = as.numeric(Copper),
      Day       = as.numeric(Day),
      Count5mm  = as.numeric(Count5mm)
    ) %>%
    filter(!is.na(Iteration), Day > 0) %>%
    mutate(Day_test = Day - sim_day_offset)
  
  # Drop still-running (incomplete) replicates: keep iterations reaching the full horizon,
  # grouped by BOTH Copper and Iteration (iteration numbers repeat per concentration).
  sim_full_day <- max(sim_raw$Day, na.rm = TRUE)
  sim_raw <- sim_raw %>%
    group_by(Copper, Iteration) %>%
    filter(max(Day) >= sim_full_day) %>%
    ungroup()
  
  # Report available concentrations and check the selected ones exist.
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
  
  used_conc <- c(control_conc, conc)
  sim_day <- sim_raw %>%
    filter(Copper %in% used_conc, Day_test %in% obs_meas_days) %>%
    group_by(Copper, Day_test) %>%
    summarise(mean_count = mean(Count5mm, na.rm = TRUE), .groups = "drop")
  
  ctrl_sim <- sim_day %>%
    filter(Copper == control_conc) %>%
    select(Day_test, ctrl = mean_count)
  
  sim_day %>%
    filter(Copper %in% conc) %>%
    left_join(ctrl_sim, by = "Day_test") %>%
    mutate(
      Treatment = exposures$Treatment[match(Copper, conc)],
      pct       = ifelse(!is.na(ctrl) & ctrl > 0, 100 * mean_count / ctrl, NA_real_),
      Source    = "Predicted"
    ) %>%
    filter(!is.na(pct)) %>%
    rename(Day = Day_test) %>%
    select(Treatment, Day, pct, Source)
}

# Full faceted figure for one basis.
build_figure <- function(src) {
  sim_pct <- build_sim_pct(src)
  
  plot_df <- bind_rows(obs_pct, sim_pct) %>%
    mutate(
      Treatment = factor(Treatment, levels = exposures$Treatment),
      Source    = factor(Source, levels = c("Observed", "Predicted"))
    )
  
  facet_labels <- setNames(
    paste0(exposures$Treatment, "  (",
           vapply(seq_len(nrow(exposures)), conc_text, character(1), src = src), ")"),
    exposures$Treatment
  )
  
  ggplot(plot_df,
         aes(x = Day, y = pct, colour = Treatment,
             linetype = Source, shape = Source,
             group = interaction(Treatment, Source))) +
    geom_hline(yintercept = 100, colour = "grey60", linewidth = 0.4, linetype = "dotted") +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.4) +
    facet_wrap(~ Treatment, ncol = 1,
               scales   = if (free_y_panels) "free_y" else "fixed",
               labeller = as_labeller(facet_labels)) +
    scale_colour_manual(values = trt_cols, guide = "none") +
    scale_linetype_manual(name = "Data", values = c(Observed = "solid", Predicted = "dashed")) +
    scale_shape_manual(name = "Data", values = c(Observed = 16, Predicted = 17)) +
    coord_cartesian(xlim = range(obs_meas_days)) +
    labs(
      title = paste0("Treatment effect over time \u2014 individuals > ", threshold_mm,
                     " mm (% of control)"),
      subtitle = paste0("Concentration basis: ",
                        if (src == "FreeIon") "free ion" else "dissolved Ni"),
      x = "Day",
      y = paste0("% of control (snails > ", threshold_mm, " mm)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      legend.position = "right",
      strip.text      = element_text(face = "bold", size = 11)
    )
}

#####################
# BUILD BOTH FIGURES (skip a basis whose sim file is missing/unreadable)

bad_src <- setdiff(ni_sources, names(sim_files))
if (length(bad_src) > 0) stop("ni_sources not in sim_files: ", paste(bad_src, collapse = ", "))

figs <- list()
for (src in ni_sources) {
  f <- tryCatch(build_figure(src),
                error = function(e) {
                  warning("Skipping ", src, ": ", conditionMessage(e)); NULL
                })
  if (!is.null(f)) figs[[src]] <- f
}
if (length(figs) == 0) stop("No figures produced — check the sim files exist.")

#####################
# OUTPUT — screen (one per basis) or PDF (one page per basis)

if (save_pdf) {
  pdf(pdf_file, width = 8, height = 10)
  for (f in figs) print(f)
  dev.off()
  cat("\nPDF saved to:", pdf_file, "  (", length(figs), "figure(s) )\n")
} else {
  for (f in figs) print(f)
}