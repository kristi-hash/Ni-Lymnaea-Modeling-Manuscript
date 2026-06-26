# TREATMENT EFFECT AS % OF CONTROL OVER TIME — cumulative reproduction
#
# Produces TWO 3-panel figures in one run: one for FreeIon, one for DissNi.
# Each figure is split into one panel per treatment (Ni1, Ni2, Ni3). For each treatment
# and timepoint, CUMULATIVE CLUTCHES are expressed as a PERCENTAGE OF THE MATCHED CONTROL:
#       value(treatment) / value(control) * 100
#   - Observed  (SOLID line + circles)   : mean treatment-tank clutches / mean Ni0-tank clutches
#   - Predicted (DASHED line + triangles): mean exposure clutches       / mean Cu=0 clutches
# 100% = same as control; < 100% reduction; > 100% increase.
#
# IMPORTANT — no Day-0 anchor: cumulative clutches start at 0 for every treatment, so the
# ratio is 0/0 (undefined) at the start. The control denominator is zero on the first
# reading, so that point is dropped; very early points (tiny control) can be noisy. The
# series therefore begin at the first day with non-zero control reproduction.
#
# Free ion and dissolved Ni are SEPARATE simulation runs; ni_sources selects the sim
# file(s) and the per-treatment concentrations matched against `Copper`. Observations are
# matched by treatment, so the observed side is identical across both figures.
#
# Metric details: observed = `Cummulative Clutches`; simulated = Cumulativeoffspring /
# eggs_per_clutch, post-reset (sim Day >= repro_reset_day; Day 180 = burn-in total).
# Time: observed anchored to first reading (repro_start_date); sim Day_test = Day - offset.
# Both series sampled ONLY at observation days (points + connecting lines); x ends at the
# last observation day. Incomplete sim replicates dropped.

rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(tibble)

#####################
# USER SETTINGS

ni_sources <- c("FreeIon", "DissNi")

sim_files <- c(
  FreeIon = "PopulationEffectlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts_FreeIon_Growth",
  DissNi  = "PopulationEffectlPredictions_Ration9500_JuvMort975_kapR47_FixedEggCosts_DissNi_Growth"
)
obs_file <- "Population_ReproductiveData.csv"

sim_day_offset   <- 180
eggs_per_clutch  <- 100
repro_reset_day  <- 181
repro_start_date <- as.Date("2023-05-26")
control_conc     <- 0
control_treatment <- "Ni0"

exposures <- tibble(
  Treatment    = c("Ni1",      "Ni2",      "Ni3"),
  free_ion     = c(0.409,      1.342,      12.114),
  dissolved_ni = c(2.30,       7.13,       65.24),
  colour       = c("#1f78b4",  "#ef8a06",  "#d7191c")
)

conc_label_metric <- "source"          # "source" | "free_ion" | "dissolved" | "both"
conc_units        <- "\u00b5g/L"

plot_day_max <- 150

# Per-treatment panels: FALSE = shared y; TRUE = independent y per panel.
free_y_panels <- FALSE

save_pdf <- FALSE
pdf_file <- "TreatmentEffect PctOfControl Repro.pdf"

#####################
# OBSERVED (computed once — matched by treatment; mean cumulative clutches per day)

obs <- read_csv(obs_file, show_col_types = FALSE)
names(obs) <- trimws(names(obs))

all_treats <- c(control_treatment, exposures$Treatment)

obs_day <- obs %>%
  filter(Treatment %in% all_treats, !is.na(`Cummulative Clutches`)) %>%
  mutate(Day = as.numeric(as.Date(Date, format = "%m/%d/%Y") - repro_start_date)) %>%
  group_by(Treatment, Day) %>%
  summarise(mean_clutch = mean(`Cummulative Clutches`), .groups = "drop")

obs_meas_days <- sort(unique(obs_day$Day))   # both series sampled at these days

ctrl_obs <- obs_day %>%
  filter(Treatment == control_treatment) %>%
  select(Day, ctrl = mean_clutch)

obs_pct <- obs_day %>%
  filter(Treatment %in% exposures$Treatment) %>%
  left_join(ctrl_obs, by = "Day") %>%
  mutate(pct = ifelse(!is.na(ctrl) & ctrl > 0, 100 * mean_clutch / ctrl, NA_real_),
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

# PREDICTED % of control for one basis (sampled at observation days, post-reset).
build_sim_pct <- function(src) {
  sim_file <- sim_files[[src]]
  conc     <- conc_for(src)
  
  sim_raw <- read_delim(sim_file, delim = " ", show_col_types = FALSE) %>%
    mutate(
      Iteration           = as.numeric(Iteration),
      Copper              = as.numeric(Copper),
      Day                 = as.numeric(Day),
      Cumulativeoffspring = as.numeric(Cumulativeoffspring)
    ) %>%
    filter(!is.na(Iteration), Day > 0) %>%
    mutate(Day_test = Day - sim_day_offset)
  
  sim_full_day <- max(sim_raw$Day, na.rm = TRUE)
  sim_raw <- sim_raw %>%
    group_by(Copper, Iteration) %>%
    filter(max(Day) >= sim_full_day) %>%
    ungroup()
  
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
    mutate(clutches = Cumulativeoffspring / eggs_per_clutch) %>%
    filter(Copper %in% used_conc, Day >= repro_reset_day, Day_test %in% obs_meas_days) %>%
    group_by(Copper, Day_test) %>%
    summarise(mean_clutch = mean(clutches, na.rm = TRUE), .groups = "drop")
  
  ctrl_sim <- sim_day %>%
    filter(Copper == control_conc) %>%
    select(Day_test, ctrl = mean_clutch)
  
  sim_day %>%
    filter(Copper %in% conc) %>%
    left_join(ctrl_sim, by = "Day_test") %>%
    mutate(
      Treatment = exposures$Treatment[match(Copper, conc)],
      pct       = ifelse(!is.na(ctrl) & ctrl > 0, 100 * mean_clutch / ctrl, NA_real_),
      Source    = "Predicted"
    ) %>%
    filter(!is.na(pct)) %>%
    rename(Day = Day_test) %>%
    select(Treatment, Day, pct, Source)
}

# Faceted figure for one basis.
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
    coord_cartesian(xlim = range(plot_df$Day)) +     # ends where observations end
    labs(
      title = "Treatment effect over time \u2014 cumulative reproduction (% of control)",
      subtitle = paste0("Concentration basis: ",
                        if (src == "FreeIon") "free ion" else "dissolved Ni"),
      x = "Day",
      y = "% of control (cumulative clutches)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(size = 10),
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
                error = function(e) { warning("Skipping ", src, ": ", conditionMessage(e)); NULL })
  if (!is.null(f)) figs[[src]] <- f
}
if (length(figs) == 0) stop("No figures produced — check the sim files exist.")

#####################
# OUTPUT

if (save_pdf) {
  pdf(pdf_file, width = 8, height = 10)
  for (f in figs) print(f)
  dev.off()
  cat("\nPDF saved to:", pdf_file, "  (", length(figs), "figure(s) )\n")
} else {
  for (f in figs) print(f)
}