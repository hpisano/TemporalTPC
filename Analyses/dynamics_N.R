# ==================== PARAMETERS =====================
library(deSolve, cubature)

# Time settings
time_start <- 0
time_end <- 100
time_step <- 0.1
times <- seq(time_start, time_end, by = time_step)

# ADD THIS: Burn-in period (in time units)
burn_in_time <- 20  # Adjust this value as needed

# Temperature parameters
temp_params <- list(
  model = "sine",
  det_params = list(mu = 20, sigma = 5, p = 24),
  noise_type = "none",
  noise_params = list(std = 1, theta = 0.1, mu = 0),
  seed = 123
)

# Growth rate parameters
growth_params <- list(
  T_opt = 25,
  epsilon = 8,
  d_inf = 0.05,
  delta_r = 1.5
)

# Population parameters
pop_params <- list(
  N0 = 10,
  alpha = 0.001
)

# ==================== COMPUTATIONS =====================

# 1. Create continuous temperature function
T_func <- do.call(generate_T, temp_params)

# 2. Create continuous growth rate function
r_func <- simple_growth_rate_r(T = T_func, 
                               T_opt = growth_params$T_opt,
                               epsilon = growth_params$epsilon,
                               d_inf = growth_params$d_inf,
                               delta_r = growth_params$delta_r)

# 3. Create continuous population systems
pop_systems <- create_population_systems(r_func, pop_params$alpha, pop_params$N0)

# Function to solve ODE and return function wrapper
solve_ode_system <- function(ode_func, y0, times, parms = list()) {
  solution <- ode(y = y0, times = times, func = ode_func, parms = parms)
  # Return a function that evaluates the ODE solution using interpolation
  # (This is the only place we need interpolation, but it's minimal)
  approxfun(solution[, 1], solution[, 2], rule = 2)
}

# Solve each model
N_dynamic_func <- solve_ode_system(pop_systems$dynamic, c(N = pop_params$N0), times)
N_null_func <- solve_ode_system(pop_systems$null, c(N = pop_params$N0), times)
N_slow_func <- solve_ode_system(pop_systems$slow, c(N = pop_params$N0), times)
K_fast_func <- pop_systems$fast  # Already a continuous function!

T_vals <- T_func(times)
r_vals <- r_func(times)
N_vals <- N_dynamic_func(times)

# Also need these for plotting:
N_null_vals <- N_null_func(times)
N_slow_vals <- N_slow_func(times)
K_fast_vals <- K_fast_func(times)
r_mean <- mean(r_vals)  

# ==================== BURN-IN PERIOD HANDLING =====================

# Create indices for burn-in period
burn_in_idx <- times <= burn_in_time
post_burn_idx <- times > burn_in_time

# Create post-burn-in time vector and values
times_post <- times[post_burn_idx]
T_vals_post <- T_vals[post_burn_idx]
r_vals_post <- r_vals[post_burn_idx]
N_vals_post <- N_vals[post_burn_idx]
N_null_vals_post <- N_null_vals[post_burn_idx]
N_slow_vals_post <- N_slow_vals[post_burn_idx]
K_fast_vals_post <- K_fast_vals[post_burn_idx]

# ==================== CONTINUOUS COMPARISON METRICS (WITH BURN-IN) =====================

abs_dev_null <- absolute_deviation(t_f = max(times),
                                   t_0 = burn_in_time,  # ADD THIS if your function supports it
                                   N_func = N_dynamic_func,
                                   m_func = N_null_func)

abs_dev_slow <- absolute_deviation(t_f = max(times),
                                   t_0 = burn_in_time,  # ADD THIS
                                   N_func = N_dynamic_func,
                                   m_func = N_slow_func)

abs_dev_fast <- absolute_deviation(t_f = max(times),
                                   t_0 = burn_in_time,  # ADD THIS
                                   N_func = N_dynamic_func,
                                   m_func = K_fast_func)

# If your absolute_deviation function doesn't support t_0 parameter, use this wrapper:
# absolute_deviation_burn <- function(t_0, t_f, N_func, m_func, ...) {
#   integrand <- function(tau) {
#     abs(N_func(tau) - m_func(tau))
#   }
#   integrate(integrand, lower = t_0, upper = t_f, ...)$value
# }

speed = abs_dev_slow / (abs_dev_fast + abs_dev_slow)

# ==================== PLOTTING (WITH BURN-IN) =====================

# Set up multi-panel plot
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

# Plot 1: Temperature over time
plot(times_post, T_vals_post, type = "l", lwd = 2, col = "red",
     xlab = "Time", ylab = "Temperature (°C)",
     main = "Temperature Dynamics (post burn-in)",
     ylim = c(min(T_vals_post) - 2, max(T_vals_post) + 2))
abline(h = temp_params$det_params$mu, lty = 3, col = "darkgray", lwd = 1.5)
abline(h = growth_params$T_opt, lty = 2, col = "blue", lwd = 1.5)
legend("topright", legend = c("Temperature", "μ (mean T)", "T_opt"), 
       col = c("red", "darkgray", "blue"), lty = c(1, 3, 2), lwd = c(2, 1.5, 1.5), bty = "n")

# Plot 2: Growth rate over time
plot(times_post, r_vals_post, type = "l", lwd = 2, col = "darkgreen",
     xlab = "Time", ylab = "Growth rate r",
     main = "Temperature-Dependent Growth Rate (post burn-in)",
     ylim = c(min(r_vals_post) - 0.1, max(r_vals_post) + 0.1))
abline(h = 0, lty = 3, col = "gray")
abline(h = mean(r_vals_post), lty = 2, col = "orange", lwd = 1.5)  # Use post-burn mean
legend("topright", legend = c("r(t)", "E[r] (post burn-in)"), 
       col = c("darkgreen", "orange"), lty = c(1, 2), lwd = c(2, 1.5), bty = "n")

# Plot 3: Population models comparison
plot(times_post, N_vals_post, type = "l", lwd = 2, col = "purple",
     xlab = "Time", ylab = "Population size N",
     main = "Population Dynamics: Model Comparison (post burn-in)",
     ylim = c(0, max(c(N_vals_post, N_null_vals_post, N_slow_vals_post, K_fast_vals_post)) * 1.1))
lines(times_post, N_null_vals_post, lwd = 2, col = "blue", lty = 2)
lines(times_post, N_slow_vals_post, lwd = 2, col = "orange", lty = 3)
lines(times_post, K_fast_vals_post, lwd = 2, col = "green", lty = 4)
abline(h = pop_params$N0, lty = 2, col = "gray", lwd = 1)
legend("bottomright", 
       legend = c("Dynamic r(t)", "Null: T=μ", "Slow: r=E[r]", "Fast: K(t)=r/α", "Initial N"),
       col = c("purple", "blue", "orange", "green", "gray"),
       lty = c(1, 2, 3, 4, 2), lwd = c(2, 2, 2, 2, 1), bty = "n")

# Add overall title
title(paste("Temperature-Dependent Population Growth: Model Comparison with Burn-in (first", burn_in_time, "time units excluded)"), outer = TRUE)

# ==================== STATISTICS AND ACF (WITH BURN-IN) =====================

# Compute CV for N (dynamic model) using post-burn-in data only
cv_N <- population_cv(t_time = times_post, N_pop = N_vals_post, weighted = TRUE)
cv_N_null <- population_cv(t_time = times_post, N_pop = N_null_vals_post, weighted = TRUE)
cv_N_slow <- population_cv(t_time = times_post, N_pop = N_slow_vals_post, weighted = TRUE)

cat("\n=== POPULATION STATISTICS COMPARISON (POST BURN-IN) ===\n")
cat("Dynamic model CV:", round(cv_N, 4), "\n")
# Compute and plot ACF for N (dynamic model) using post-burn-in data
acf_N <- population_acf(t_time = times_post, N_pop = N_vals_post, 
                        tau_max = 25, n_lags = 100)

# Print results
cat("\n=== DEVIATION METRICS (WITH BURN-IN) ===\n")
cat("Absolute deviation from null model:", round(abs_dev_null, 4), "\n")
cat("Absolute deviation from slow model:", round(abs_dev_slow, 4), "\n")
cat("Absolute deviation from fast model:", round(abs_dev_fast, 4), "\n")
cat("Speed metric:", round(speed, 4), "\n")