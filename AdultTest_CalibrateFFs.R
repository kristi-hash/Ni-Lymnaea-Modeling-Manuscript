# =============================================================================
# Lymnaea stagnalis DEB-IBM — Adult food factor calibration
# Two-phase calibration: FF1 (rearing, before FF switch) and FF2 (exposure)
# Targets: PopulationWetWeight at Day 151 and Day 179
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
sim_file  <- "DynamicEggCost_cvChange_CalibrateAdultFF2_FFChange137_V2"
n_best    <- 100

obs_WW151 <- 2.0237   # g, observed mean wet weight at Day 151
obs_WW179 <- 2.2901   # g, observed mean wet weight at Day 179

# Days used to extract food factors
# FF1: well within rearing phase (before switch)
# FF2: well within exposure phase (after switch)
day_FF1 <- 100
day_FF2 <- 160

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------
df <- read.table(sim_file, header = TRUE, sep = " ", fill = TRUE) %>%
  filter(Iteration != "Iteration") %>%
  mutate(across(c(Iteration, Day, FoodFactor, PopulationWetWeight), as.numeric)) %>%
  filter(!is.na(Iteration), Iteration >= 1, !is.na(Day))

cat("Iterations:", n_distinct(df$Iteration), "\n")
cat("Day range:", min(df$Day), "-", max(df$Day), "\n\n")

# -----------------------------------------------------------------------------
# 3. Extract FF1, FF2, WW151, WW179 per iteration
# -----------------------------------------------------------------------------
ff1   <- df %>% filter(Day == day_FF1) %>%
  group_by(Iteration) %>% slice_max(FoodFactor, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(Iteration, FF1 = FoodFactor)

ff2   <- df %>% filter(Day == day_FF2) %>%
  group_by(Iteration) %>% slice_max(FoodFactor, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(Iteration, FF2 = FoodFactor)

ww151 <- df %>% filter(Day == 151) %>%
  group_by(Iteration) %>% slice_max(PopulationWetWeight, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(Iteration, WW151 = PopulationWetWeight)

ww179 <- df %>% filter(Day == 179) %>%
  group_by(Iteration) %>% slice_max(PopulationWetWeight, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(Iteration, WW179 = PopulationWetWeight)

dat <- ff1 %>%
  inner_join(ff2,   by = "Iteration") %>%
  inner_join(ww151, by = "Iteration") %>%
  inner_join(ww179, by = "Iteration") %>%
  filter(WW151 > 0, WW179 > 0)

cat("Valid iterations:", nrow(dat), "\n")
cat(sprintf("FF1 range: %.4f - %.4f\n", min(dat$FF1), max(dat$FF1)))
cat(sprintf("FF2 range: %.4f - %.4f\n", min(dat$FF2), max(dat$FF2)))
cat(sprintf("WW151 range: %.4f - %.4f g\n", min(dat$WW151), max(dat$WW151)))
cat(sprintf("WW179 range: %.4f - %.4f g\n\n", min(dat$WW179), max(dat$WW179)))

# -----------------------------------------------------------------------------
# 4a. Joint calibration: SSE on both WW151 and WW179
# -----------------------------------------------------------------------------
dat <- dat %>%
  mutate(SSE_joint = (WW151 - obs_WW151)^2 + (WW179 - obs_WW179)^2,
         SSE_151   = (WW151 - obs_WW151)^2)

joint   <- dat %>% arrange(SSE_joint)
best_j  <- joint %>% slice_head(n = 1)
top_j   <- joint %>% slice_head(n = n_best)

cat("=== Joint calibration (WW151 + WW179) ===\n")
cat("Best single pair:\n")
cat(sprintf("  FF1:   %.4f\n", best_j$FF1))
cat(sprintf("  FF2:   %.4f\n", best_j$FF2))
cat(sprintf("  WW151: %.4f g  (obs: %.4f g, diff: %+.4f g)\n",
            best_j$WW151, obs_WW151, best_j$WW151 - obs_WW151))
cat(sprintf("  WW179: %.4f g  (obs: %.4f g, diff: %+.4f g)\n",
            best_j$WW179, obs_WW179, best_j$WW179 - obs_WW179))
cat(sprintf("  SSE:   %.6f\n\n", best_j$SSE_joint))

cat(sprintf("Median of top %d:\n", n_best))
cat(sprintf("  FF1: median=%.4f, 95%%CI=[%.4f, %.4f]\n",
            median(top_j$FF1), quantile(top_j$FF1, 0.025), quantile(top_j$FF1, 0.975)))
cat(sprintf("  FF2: median=%.4f, 95%%CI=[%.4f, %.4f]\n",
            median(top_j$FF2), quantile(top_j$FF2, 0.025), quantile(top_j$FF2, 0.975)))
cat(sprintf("  WW151: mean=%.4f g, sd=%.4f g  (obs: %.4f g)\n",
            mean(top_j$WW151), sd(top_j$WW151), obs_WW151))
cat(sprintf("  WW179: mean=%.4f g, sd=%.4f g  (obs: %.4f g)\n\n",
            mean(top_j$WW179), sd(top_j$WW179), obs_WW179))

n_both <- dat %>%
  filter(between(WW151, obs_WW151 - 0.05, obs_WW151 + 0.05),
         between(WW179, obs_WW179 - 0.05, obs_WW179 + 0.05)) %>%
  nrow()
cat(sprintf("Iterations matching both targets (±0.05 g): %d\n\n", n_both))

# -----------------------------------------------------------------------------
# 4b. FF1-only calibration: SSE on WW151 only
# -----------------------------------------------------------------------------
ff1_cal  <- dat %>% arrange(SSE_151)
best_f1  <- ff1_cal %>% slice_head(n = 1)
top_f1   <- ff1_cal %>% slice_head(n = n_best)

cat("=== FF1 calibration (WW151 only) ===\n")
cat("Best single FF1:\n")
cat(sprintf("  FF1:   %.4f  (FF2=%.4f — not used in this SSE)\n",
            best_f1$FF1, best_f1$FF2))
cat(sprintf("  WW151: %.4f g  (obs: %.4f g, diff: %+.4f g)\n\n",
            best_f1$WW151, obs_WW151, best_f1$WW151 - obs_WW151))

cat(sprintf("Median of top %d FF1:\n", n_best))
cat(sprintf("  FF1: median=%.4f, 95%%CI=[%.4f, %.4f]\n",
            median(top_f1$FF1), quantile(top_f1$FF1, 0.025), quantile(top_f1$FF1, 0.975)))
cat(sprintf("  WW151: mean=%.4f g, sd=%.4f g  (obs: %.4f g)\n\n",
            mean(top_f1$WW151), sd(top_f1$WW151), obs_WW151))

# -----------------------------------------------------------------------------
# 5. Plots
# -----------------------------------------------------------------------------
p1 <- ggplot(dat, aes(x = FF1, y = WW151)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "steelblue") +
  geom_point(data = top_j, colour = "darkorange", size = 1.5, alpha = 0.8) +
  geom_point(data = best_j, colour = "red", size = 3, shape = 8) +
  geom_hline(yintercept = obs_WW151, linetype = "dashed", colour = "red") +
  labs(x = "FF1 (rearing)", y = "Predicted WW Day 151 (g)",
       caption = "Orange = joint top 100; red star = best; dashed = observed") +
  theme_bw(base_size = 11)

p2 <- ggplot(dat, aes(x = FF2, y = WW179)) +
  geom_point(alpha = 0.15, size = 0.8, colour = "steelblue") +
  geom_point(data = top_j, colour = "darkorange", size = 1.5, alpha = 0.8) +
  geom_point(data = best_j, colour = "red", size = 3, shape = 8) +
  geom_hline(yintercept = obs_WW179, linetype = "dashed", colour = "red") +
  labs(x = "FF2 (exposure)", y = "Predicted WW Day 179 (g)",
       caption = "Orange = joint top 100; red star = best; dashed = observed") +
  theme_bw(base_size = 11)

p3 <- ggplot(top_j, aes(x = WW151, y = WW179, colour = SSE_joint)) +
  geom_point(size = 2.5) +
  geom_point(data = best_j, colour = "red", size = 3, shape = 8) +
  geom_vline(xintercept = obs_WW151, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = obs_WW179, linetype = "dashed", colour = "red") +
  scale_colour_viridis_c(direction = -1) +
  labs(x = "Predicted WW Day 151 (g)", y = "Predicted WW Day 179 (g)",
       title = sprintf("Top %d: WW151 vs WW179", n_best)) +
  theme_bw(base_size = 11)

p4 <- ggplot(top_j, aes(x = FF1, y = FF2, colour = SSE_joint)) +
  geom_point(size = 2.5) +
  geom_point(data = best_j, colour = "red", size = 3, shape = 8) +
  scale_colour_viridis_c(direction = -1) +
  labs(x = "FF1 (rearing)", y = "FF2 (exposure)",
       title = sprintf("Top %d: FF1 vs FF2", n_best)) +
  theme_bw(base_size = 11)

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

cat(sprintf("Top 20 joint pairs:\n"))
print(joint %>% slice_head(n = 20) %>%
        select(Iteration, FF1, FF2, WW151, WW179, SSE_joint) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))),
      row.names = FALSE)