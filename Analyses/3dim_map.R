# Load required libraries
library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)
library(future.apply)
library(parallel)

# ==================== USER-DEFINED PARAMETERS =====================
# Set these values at the beginning - they will be used throughout

# Time settings
t_prime_start <- 0
t_prime_end <- 26
t_prime_step <- 0.01
t_prime <- seq(t_prime_start, t_prime_end, by = t_prime_step)
burn_in_time <- 25

# z (normalized temperature) parameters
z_params <- list(
  model = "sine",    # Change if needed
  phase = 0          # Change if needed
)

# u (thermal response) parameters
d_inf <- 0.1         # Set d_inf value here
# Note: P_offset, P_amp will be varied in the loop

# Normalized population parameters
N_prime0 <- 0.5      # Set initial population value here
# Note: P_time will be varied in the loop

# Parameter ranges for the sweep
P_amp_vals <- seq(0, 2, length.out = 5)      
P_offset_vals <- seq(-4, 1, length.out = 8)  
P_time_vals <- 10^seq(-2, 2, length.out = 5)     

# ==================== PARAMETER SWEEP =====================
# Run simulations for all parameter combinations in parallel

total_combinations <- length(P_offset_vals) * length(P_amp_vals) * length(P_time_vals)
cat("Running", total_combinations, "simulations in parallel...\n")
cat("Parameters:\n")
cat("  d_inf =", d_inf, "\n")
cat("  N_prime0 =", N_prime0, "\n")
cat("  z model =", z_params$model, "\n")
cat("  z phase =", z_params$phase, "\n\n")

# Set up parallel processing
plan(multisession, workers = 1)

# Create all parameter combinations
param_combinations <- expand.grid(
  P_offset = P_offset_vals,
  P_amp = P_amp_vals,
  P_time = P_time_vals,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Run simulations in parallel
results_list <- future_lapply(seq_len(nrow(param_combinations)), function(i) {
  run_range_3dim_simulation(
    param_combinations$P_offset[i],
    param_combinations$P_amp[i],
    param_combinations$P_time[i]
  )
}, future.seed = TRUE)



# Combine all results
results_df <- do.call(rbind, results_list)
# ==================== CREATE HEATMAPS =====================

# Create heatmaps for E[N_prime] (10 subplots, one for each P_offset)
E_N_prime_plots <- list()
for (i in seq_along(P_offset_vals)) {
  offset_val <- P_offset_vals[i]
  subset_data <- results_df[results_df$P_offset == offset_val, ]
  
  p <- create_heatmap(
    data = subset_data,
    metric = "E_N_prime",
    title_prefix = "E[N'] after burn-in",
    color_palette = "plasma"
  )
  
  E_N_prime_plots[[i]] <- p
}

# Create heatmaps for speed 
speed_plots <- list()
for (i in seq_along(P_offset_vals)) {
  offset_val <- P_offset_vals[i]
  subset_data <- results_df[results_df$P_offset == offset_val, ]
  
  p <- create_heatmap(
    data = subset_data,
    metric = "speed",
    title_prefix = "Speed metric",
    color_palette = "magma"
  )
  
  speed_plots[[i]] <- p
}
# ==================== DISPLAY RESULTS =====================

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  E_N_prime_facet <- ggplot(results_df, aes(x = P_amp, y = P_time, fill = E_N_prime)) +
    geom_tile() +
    scale_fill_viridis(option = "plasma", limits = c(0, 1)) +
    facet_wrap(~ round(P_offset, 1), nrow = 2) +
    labs(title = paste("E[N'] after burn-in\n(d_inf =", d_inf, 
                       ", N_prime0 =", N_prime0, ")"),
         x = "P_amp", y = "P_time", fill = "E[N']") +
    theme_minimal()
  
  # Create a new column for conditional coloring
  results_df$speed_display <- ifelse(
    is.na(results_df$E_N_prime) | results_df$E_N_prime == 0,
    NA,  # Will be colored gray via na.value
    results_df$speed
  )
  
  # Modified speed plot to show gray squares for missing/zero E[N']
  speed_facet <- ggplot(results_df, aes(x = P_amp, y = P_time, fill = speed_display)) +
    geom_tile() +
    scale_fill_viridis(
      option = "magma", 
      limits = c(0, 1), 
      na.value = "gray70",  # This will color NA values gray
      guide = guide_colorbar(title = "Speed")
    ) +
    facet_wrap(~ round(P_offset, 1), nrow = 2) +
    labs(title = paste("Speed metric\n(d_inf =", d_inf, 
                       ", N_prime0 =", N_prime0, ")"),
         x = "P_amp", y = "P_time") +
    theme_minimal()
  
  # Display faceted plots
  print(E_N_prime_facet)
  print(speed_facet)
}


# ==================== CREATE ABS_DEV_SLOW HEATMAP (LOG SCALE) =====================
cat("\n=== Heatmap of abs_dev_slow (log scale) ===\n")

# Prepare data for plotting - using log scale for abs_dev_slow
results_df$abs_dev_slow_log <- log10(results_df$abs_dev_slow + 1)  # Add 1 to avoid log(0)
results_df$abs_dev_slow_log_display <- ifelse(
  is.na(results_df$E_N_prime) | results_df$E_N_prime == 0,
  NA,
  results_df$abs_dev_slow_log
)

# Calculate reasonable limits for the color scale
abs_dev_slow_min <- min(results_df$abs_dev_slow_log_display, na.rm = TRUE)
abs_dev_slow_max <- max(results_df$abs_dev_slow_log_display, na.rm = TRUE)

abs_dev_slow_plot <- ggplot(results_df, aes(x = P_amp, y = P_time, fill = abs_dev_slow_log_display)) +
  geom_tile() +
  scale_fill_viridis(
    option = "viridis", 
    na.value = "gray70",
    name = "log10(abs_dev_slow + 1)",
    limits = c(abs_dev_slow_min, abs_dev_slow_max)  # Use calculated limits
  ) +
  facet_wrap(~ round(P_offset, 1), nrow = 2) +
  labs(
    title = paste("Absolute deviation from slow system (log scale)\n(d_inf =", d_inf, 
                  ", N_prime0 =", N_prime0, ")"),
    x = "P_amp", 
    y = "P_time"
  ) +
  theme_minimal()

print(abs_dev_slow_plot)

# ==================== CREATE ABS_DEV_FAST HEATMAP (LOG SCALE + 1) =====================
cat("\n=== Heatmap of abs_dev_fast (log10(value + 1) scale) ===\n")

# Prepare data for plotting - using log10(value + 1) scale for abs_dev_fast
results_df$abs_dev_fast_log <- log10(results_df$abs_dev_fast + 1)  # Add 1 as requested
results_df$abs_dev_fast_log_display <- ifelse(
  is.na(results_df$E_N_prime) | results_df$E_N_prime == 0,
  NA,
  results_df$abs_dev_fast_log
)

# Calculate reasonable limits for the color scale
abs_dev_fast_min <- min(results_df$abs_dev_fast_log_display, na.rm = TRUE)
abs_dev_fast_max <- max(results_df$abs_dev_fast_log_display, na.rm = TRUE)

abs_dev_fast_plot <- ggplot(results_df, aes(x = P_amp, y = P_time, fill = abs_dev_fast_log_display)) +
  geom_tile() +
  scale_fill_viridis(
    option = "plasma", 
    na.value = "gray70",
    name = "log10(abs_dev_fast + 1)",
    limits = c(abs_dev_fast_min, abs_dev_fast_max)  # Use calculated limits
  ) +
  facet_wrap(~ round(P_offset, 1), nrow = 2) +
  labs(
    title = paste("Absolute deviation from fast system (log10(value + 1) scale)\n(d_inf =", d_inf, 
                  ", N_prime0 =", N_prime0, ")"),
    x = "P_amp", 
    y = "P_time"
  ) +
  theme_minimal()

print(abs_dev_fast_plot)

# ==================== CREATE MATCH MAP =====================
cat("\n=== Match Map (which deviation is bigger?) ===\n")

# Prepare data for plotting
results_df$dev_comparison <- ifelse(
  is.na(results_df$E_N_prime) | results_df$E_N_prime == 0,
  NA,
  ifelse(results_df$abs_dev_slow > results_df$abs_dev_fast,
         "slow > fast",
         "fast ≥ slow")
)

match_map_plot <- ggplot(results_df, aes(x = P_amp, y = P_time, fill = dev_comparison)) +
  geom_tile() +
  scale_fill_manual(
    values = c("slow > fast" = "blue", 
               "fast ≥ slow" = "red"), 
    na.value = "gray70",
    name = "Which is bigger?"
  ) +
  facet_wrap(~ round(P_offset, 1), nrow = 2) +
  labs(
    title = paste("Match Map: Which deviation dominates?\n(d_inf =", d_inf, 
                  ", N_prime0 =", N_prime0, ")"),
    x = "P_amp", 
    y = "P_time"
  ) +
  theme_minimal()

print(match_map_plot)