# FOOD FACTOR CALIBRATION (CORRECT)
# Version 2 has updated output processing (writes away 95% confidence intervals of posterior distribution)
# Version 2 also includes density plots of the parameters

# Clear Memory
rm(list = ls())

# Load packages
library(readr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(cowplot)

#####################
# OPTIMIZE FOOD FACTOR

# Read simulation output
#sim <- read_delim("FoodCorrV2_TKV6_Mattson_CalElimination_18Degrees", delim = " ")
sim <- read_delim("FoodCorrV2_TKV6_Mattson_CalElimination_26Degrees", delim = " ")

# Keep only Day 45
#sim <- sim[which(sim$Day == 48), ]
sim <- sim[which(sim$Day == 52), ]
#sim <- sim[which(sim$Copper == 0), ]

# Convert InternalConcentration to numeric
sim$InternalConcentration <- as.numeric(sim$InternalConcentration)

# Calculate SSE for growth relative to target weight of 0.35g
#sim$SSE_growth <- (sim$InternalConcentration - 0.47066778)^2 # 18 Degrees
sim$SSE_growth <- (sim$InternalConcentration - 0.118666936)^2 # 26 Degrees

# Create unique ID for each iteration/parameter combination
sim$SSEid <- paste(sim$Iteration, sim$EliminationRate, sim$Copper)

# Summed SSE across individuals for each parameter set
SSE_sum_growth <- tapply(X = sim$SSE_growth, INDEX = sim$SSEid, sum)

# Compute log-likelihood (simplified)
LL_sum <- (-1/2) * log(SSE_sum_growth)

# Sort by log-likelihood (best-fitting parameter sets first)
LL_sum_order <- sort(LL_sum, decreasing = TRUE)

# Select 100 best fits and extract parameter values
best <- head(LL_sum_order, n = 100)
selection <- strsplit(rownames(best), " ")
df <- data.frame(matrix(unlist(selection), nrow = 100, byrow = TRUE))

# Clean and convert EliminationRate (assumed to be column 2)
df$EliminationRate <- as.numeric(as.character(df$X2))

# Compute summary statistics
Assim_cor_max     <- max(df$EliminationRate, na.rm = TRUE)
Assim_cor_median  <- median(df$EliminationRate, na.rm = TRUE)
Assim_cor_Q1      <- quantile(df$EliminationRate, 0.025, na.rm = TRUE)
Assim_cor_Q3      <- quantile(df$EliminationRate, 0.975, na.rm = TRUE)
LL_mean           <- mean(best)
LL_max            <- max(best)

# Print summary
cat("Summary of Food Factor from Top 100 Fits:\n")
cat("  Max:     ", Assim_cor_max, "\n")
cat("  Median:  ", Assim_cor_median, "\n")
cat("  2.5% CI: ", Assim_cor_Q1, "\n")
cat("  97.5% CI:", Assim_cor_Q3, "\n")
cat("  Mean LL: ", LL_mean, "\n")
cat("  Max LL:  ", LL_max, "\n")

# Density plot of EliminationRate for best 100
ggplot(df, aes(x = EliminationRate)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = Assim_cor_median, linetype = "dashed", color = "darkblue") +
  geom_vline(xintercept = Assim_cor_Q1, linetype = "dotted", color = "red") +
  geom_vline(xintercept = Assim_cor_Q3, linetype = "dotted", color = "red") +
  labs(
    title = "Posterior Density of Food Factor (Top 100 Fits)",
    x = "Food Factor",
    y = "Density"
  ) +
  theme_minimal()
