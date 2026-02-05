# Load required libraries


library(deSolve)
library(ggplot2)
library(gridExtra)
library(viridis)

# ==================== USER-DEFINED PARAMETERS =====================
# Set these values at the beginning - they will be used throughout

# Time settings
t_prime_start <- 0
t_prime_end <- 20
t_prime_step <- 0.01
t_prime <- seq(t_prime_start, t_prime_end, by = t_prime_step)
burn_in_time <- 10

# z (normalized temperature) parameters
z_params <- list(
  model = "sine",    # Change if needed
  phase = 0          # Change if needed
)

# u (thermal response) parameters
d_inf <- 0.2         # Set d_inf value here
# Note: P_offset, P_amp will be varied in the loop

# Normalized population parameters
N_prime0 <- 0.5      # Set initial population value here
# Note: P_time will be varied in the loop

# Parameter ranges for the sweep
P_amp_vals <- seq(0, 1.5, length.out = 5)      
P_offset_vals <- seq(-2.5, 0.5  , length.out = 8)  
P_time_vals <- seq(1, 40, length.out = 5)     




# ==================== PARAMETER SWEEP =====================
# Run simulations for all parameter combinations
results <- list()
total_combinations <- length(P_offset_vals) * length(P_amp_vals) * length(P_time_vals)
cat("Running", total_combinations, "simulations...\n")
cat("Parameters:\n")
cat("  d_inf =", d_inf, "\n")
cat("  N_prime0 =", N_prime0, "\n")
cat("  z model =", z_params$model, "\n")
cat("  z phase =", z_params$phase, "\n\n")

counter <- 0
for (i in seq_along(P_offset_vals)) {
  for (j in seq_along(P_amp_vals)) {
    for (k in seq_along(P_time_vals)) {
      counter <- counter + 1
      if (counter %% 100 == 0) {
        cat("  Completed", counter, "of", total_combinations, "simulations\n")
      }
      
      result <- run_range_3dim_simulation(P_offset_vals[i], P_amp_vals[j], P_time_vals[k])
      results[[counter]] <- result
    }
  }
}

# Combine all results
results_df <- do.call(rbind, results)

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
  
  speed_facet <- ggplot(results_df, aes(x = P_amp, y = P_time, fill = speed)) +
    geom_tile() +
    scale_fill_viridis(option = "magma", limits = c(0, 1)) +
    facet_wrap(~ round(P_offset, 1), nrow = 2) +
    labs(title = paste("Speed metric\n(d_inf =", d_inf, 
                       ", N_prime0 =", N_prime0, ")"),
         x = "P_amp", y = "P_time", fill = "Speed") +
    theme_minimal()
  
  # Display faceted plots
  print(E_N_prime_facet)
  print(speed_facet)
}
