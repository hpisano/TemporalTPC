# ==================== PARAMETERS =====================

library(deSolve, cubature)

# Normalized time settings
t_prime_start <- 0
t_prime_end <- 10
t_prime_step <- 0.01
t_prime <- seq(t_prime_start, t_prime_end, by = t_prime_step)

# z (normalized temperature) parameters
z_params <- list(
  model = "sine",
  period = 1,  # Period in normalized time units
  phase = 0
)

# u (thermal response) parameters
u_params <- list(
  P_offset = 0.5,
  P_amp = 1.2,
  d_inf = 0.1
)

# Normalized population parameters
N_prime_params <- list(
  N_prime0 = 0.1,
  P_time = 1.5
)

# ==================== COMPUTATIONS =====================

# 1. Create continuous z function (normalized temperature)
z_func <- create_z_function(model = z_params$model,
                            period = z_params$period,
                            phase = z_params$phase)

# 2. Create continuous u function (thermal response)
u_func <- growth_rate_3dim_u(P_offset = u_params$P_offset,
                             P_amp = u_params$P_amp,
                             d_inf = u_params$d_inf,
                             z = z_func)

# 3. Create continuous population systems (normalized)
pop_systems_prime <- create_normalized_systems(u_func, N_prime_params$P_time, N_prime_params$N_prime0)

# 4. Solve all ODE systems
library(deSolve)

# Function to solve ODE and return function wrapper
solve_ode_system <- function(ode_func, y0, times, parms = list()) {
  solution <- ode(y = y0, times = times, func = ode_func, parms = parms)
  approxfun(solution[, 1], solution[, 2], rule = 2)
}

# Solve each model
N_prime_dynamic_func <- solve_ode_system(pop_systems_prime$dynamic, c(N_prime = N_prime_params$N_prime0), t_prime)
N_prime_null_func <- solve_ode_system(pop_systems_prime$null, c(N_prime = N_prime_params$N_prime0), t_prime)
N_prime_slow_func <- solve_ode_system(pop_systems_prime$slow, c(N_prime = N_prime_params$N_prime0), t_prime)
K_prime_fast_func <- pop_systems_prime$fast  # Already a continuous function!

# Evaluate functions at time points for plotting
z_vals <- z_func(t_prime)
u_vals <- u_func(t_prime)
N_prime_vals <- N_prime_dynamic_func(t_prime)
N_prime_null_vals <- N_prime_null_func(t_prime)
N_prime_slow_vals <- N_prime_slow_func(t_prime)
K_prime_fast_vals <- K_prime_fast_func(t_prime)
u_mean <- mean(u_vals)  # For reference line in plot

# ==================== CONTINUOUS COMPARISON METRICS =====================

abs_dev_null_prime <- absolute_deviation(t_f = max(t_prime),
                                        N_func = N_prime_dynamic_func,
                                        m_func = N_prime_null_func)

abs_dev_slow_prime <- absolute_deviation(t_f = max(t_prime),
                                        N_func = N_prime_dynamic_func,
                                        m_func = N_prime_slow_func)

abs_dev_fast_prime <- absolute_deviation(t_f = max(t_prime),
                                        N_func = N_prime_dynamic_func,
                                        m_func = K_prime_fast_func)

# ==================== PLOTTING =====================

# Set up multi-panel plot
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

# Plot 1: Normalized temperature (z) over normalized time
plot(t_prime, z_vals, type = "l", lwd = 2, col = "blue",
     xlab = "Normalized time (t')", ylab = "z (normalized temperature)",
     main = "Normalized Temperature Dynamics",
     ylim = c(min(z_vals) - 0.5, max(z_vals) + 0.5))
abline(h = 0, lty = 3, col = "darkgray", lwd = 1.5)
legend("topright", legend = c("z(t')", "zero (μ_z)"), 
       col = c("blue", "darkgray"), lty = c(1, 3), lwd = c(2, 1.5), bty = "n")

# Plot 2: Thermal response (u) over normalized time
plot(t_prime, u_vals, type = "l", lwd = 2, col = "darkorange",
     xlab = "Normalized time (t')", ylab = "u (thermal response)",
     main = "Thermal Response Function",
     ylim = c(min(u_vals) - 0.1, max(u_vals) + 0.1))
abline(h = 0, lty = 3, col = "gray")
abline(h = u_mean, lty = 2, col = "red", lwd = 1.5)
legend("topright", legend = c("u(t')", "E[u]"), 
       col = c("darkorange", "red"), lty = c(1, 2), lwd = c(2, 1.5), bty = "n")

# Plot 3: Normalized population models comparison
plot(t_prime, N_prime_vals, type = "l", lwd = 2, col = "purple",
     xlab = "Normalized time (t')", ylab = "N' (normalized population)",
     main = "Normalized Population Dynamics: Model Comparison",
     ylim = c(0, max(c(N_prime_vals, N_prime_null_vals, N_prime_slow_vals, K_prime_fast_vals)) * 1.1))
lines(t_prime, N_prime_null_vals, lwd = 2, col = "blue", lty = 2)
lines(t_prime, N_prime_slow_vals, lwd = 2, col = "orange", lty = 3)
lines(t_prime, K_prime_fast_vals, lwd = 2, col = "green", lty = 4)
abline(h = N_prime_params$N_prime0, lty = 2, col = "gray", lwd = 1)
abline(h = 1, lty = 3, col = "red", lwd = 1.5)
legend("right", 
       legend = c("Dynamic u(t')", "Null: z=0", "Slow: u=E[u]", "Fast: u/P_time", "Initial N'", "Carrying cap."),
       col = c("purple", "blue", "orange", "green", "gray", "red"),
       lty = c(1, 2, 3, 4, 2, 3), lwd = c(2, 2, 2, 2, 1, 1.5), bty = "n")

# Add overall title
title("Normalized Model: z, u, and N' with Model Comparison", outer = TRUE)

# ==================== STATISTICS AND ACF =====================

# Compute CV for N' (dynamic model)
cv_Nprime <- population_cv(t_time = t_prime, N_pop = N_prime_vals, weighted = TRUE)
cv_Nprime_null <- population_cv(t_time = t_prime, N_pop = N_prime_null_vals, weighted = TRUE)
cv_Nprime_slow <- population_cv(t_time = t_prime, N_pop = N_prime_slow_vals, weighted = TRUE)

cat("\n=== NORMALIZED POPULATION STATISTICS COMPARISON ===\n")
cat("Dynamic model CV:", round(cv_Nprime, 4), "\n")

# Compute and plot ACF for N' (dynamic model)
acf_Nprime <- population_acf(t_time = t_prime, N_pop = N_prime_vals, 
                             tau_max = 2.5, n_lags = 100)

# Print results
cat("\n=== CONTINUOUS MODEL COMPARISON (NORMALIZED) ===\n")
cat("Null model (z constant at 0):\n")
cat("  Absolute Deviation:", round(abs_dev_null_prime, 4), "\n")

cat("Slow model (u constant at mean):\n")
cat("  Absolute Deviation:", round(abs_dev_slow_prime, 4), "\n")

cat("Fast model (u/P_time):\n")
cat("  Absolute Deviation:", round(abs_dev_fast_prime, 4), "\n")