# ==================== PARAMETERS =====================
library(deSolve, cubature)

# Normalized time settings
t_prime_start <- 0
t_prime_end <- 20
t_prime_step <- 0.01
t_prime <- seq(t_prime_start, t_prime_end, by = t_prime_step)

# Burn-in period (in normalized time units)
burn_in_time <- 10

# z (normalized temperature) parameters
z_params <- list(
  model = "sine",
  phase = 0
)

# u (thermal response) parameters
u_params <- list(
  P_offset = 0.5,
  P_amp = 2,
  d_inf = 0.1
)

# Normalized population parameters
N_prime_params <- list(
  N_prime0 = 0.5,
  P_time = 10
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

# Function to solve ODE and return function wrapper
solve_ode_system <- function(ode_func, y0, times, parms = list()) {
  solution <- ode(y = y0, times = times, func = ode_func, parms = parms)
  approxfun(solution[, 1], solution[, 2], rule = 2)
}

# Solve each model
N_prime_dynamic_func <- solve_ode_system(pop_systems_prime$dynamic, c(N_prime = N_prime_params$N_prime0), t_prime)
N_prime_null_func <- solve_ode_system(pop_systems_prime$null, c(N_prime = N_prime_params$N_prime0), t_prime)
N_prime_slow_func <- solve_ode_system(pop_systems_prime$slow, c(N_prime = N_prime_params$N_prime0), t_prime)
K_prime_fast_func <- pop_systems_prime$fast

# Evaluate functions at time points for plotting
z_vals <- z_func(t_prime)
u_vals <- u_func(t_prime)
N_prime_vals <- N_prime_dynamic_func(t_prime)
N_prime_null_vals <- N_prime_null_func(t_prime)
N_prime_slow_vals <- N_prime_slow_func(t_prime)
K_prime_fast_vals <- K_prime_fast_func(t_prime)

# Apply threshold to prevent extremely small population values
N_prime_vals[N_prime_vals < 0.01] <- 0
N_prime_null_vals[N_prime_null_vals < 0.01] <- 0
N_prime_slow_vals[N_prime_slow_vals < 0.01] <- 0

# ==================== BURN-IN PERIOD HANDLING =====================

# Create indices for burn-in period
burn_in_idx <- t_prime <= burn_in_time
post_burn_idx <- t_prime > burn_in_time

# Create post-burn-in time vector and values
t_prime_post <- t_prime[post_burn_idx]
z_vals_post <- z_vals[post_burn_idx]
u_vals_post <- u_vals[post_burn_idx]
N_prime_vals_post <- N_prime_vals[post_burn_idx]
N_prime_null_vals_post <- N_prime_null_vals[post_burn_idx]
N_prime_slow_vals_post <- N_prime_slow_vals[post_burn_idx]
K_prime_fast_vals_post <- K_prime_fast_vals[post_burn_idx]

# ==================== CONTINUOUS COMPARISON METRICS (WITH BURN-IN) =====================

# Create functions that only exist after burn-in period for metrics
abs_dev_null_prime <- absolute_deviation(t_f = max(t_prime),
                                        t_0 = burn_in_time,
                                        N_func = N_prime_dynamic_func,
                                        m_func = N_prime_null_func)

abs_dev_slow_prime <- absolute_deviation(t_f = max(t_prime),
                                        t_0 = burn_in_time,
                                        N_func = N_prime_dynamic_func,
                                        m_func = N_prime_slow_func)

abs_dev_fast_prime <- absolute_deviation(t_f = max(t_prime),
                                        t_0 = burn_in_time,
                                        N_func = N_prime_dynamic_func,
                                        m_func = K_prime_fast_func)

speed = abs_dev_slow_prime / (abs_dev_fast_prime + abs_dev_slow_prime)
# ==================== PLOTTING (WITH BURN-IN AND PARAMETERS IN TITLE) =====================

# Create parameter string for title
param_string <- paste0(
  "  P_offset=", u_params$P_offset, 
  ",  P_amp=", u_params$P_amp, 
  ",  P_time=", N_prime_params$P_time, 
  ",  d_inf=", u_params$d_inf, 
  ",  burn-in=", burn_in_time
)

# Create performance metrics string for subtitle
metrics_string <- paste0(
  "null_dev=", round(abs_dev_null_prime, 3), 
  ", slow_dev=", round(abs_dev_slow_prime, 3), 
  ", fast_dev=", round(abs_dev_fast_prime, 3), 
  ", speed=", round(speed, 3)
)

# Set up multi-panel plot with appropriate margins
par(mfrow = c(3, 1), mar = c(4, 4, 3, 10), oma = c(2, 0, 6, 0))

# Plot 1: Normalized temperature (z) over normalized time
plot(t_prime_post, z_vals_post, type = "l", lwd = 2, col = "#cb1f28",
     xlab = "Normalized time (t')", ylab = "z (normalized temperature)",
     main = "Normalized Temperature Dynamics (post burn-in)",
     ylim = c(min(z_vals_post) - 0.5, max(z_vals_post) + 0.5))
abline(h = 0, lty = 3, col = "darkgray", lwd = 1.5)
# Legend positioned to the right side but within the plot area
legend(x = max(t_prime_post) * 1.02, y = max(z_vals_post) * 0.9,
       legend = c("z(t')", "zero (μ_z)"), 
       col = c("#cb1f28", "darkgray"), lty = c(1, 3), lwd = c(2, 1.5), 
       bty = "n", xpd = TRUE, cex = 0.9)

# Plot 2: Thermal response (u) over normalized time
plot(t_prime_post, u_vals_post, type = "l", lwd = 2, col = "#4e6fc1",
     xlab = "Normalized time (t')", ylab = "u (thermal response)",
     main = "Thermal Response Function (post burn-in)",
     ylim = c(min(u_vals_post) - 0.1, max(u_vals_post) + 0.1))
abline(h = 0, lty = 3, col = "gray")
abline(h = mean(u_vals_post), lty = 3, col = "#cd5cbf", lwd = 1.5)
# Legend positioned to the right side but within the plot area
legend(x = max(t_prime_post) * 1.02, y = max(u_vals_post) * 0.9,
       legend = c("u(t')", "E[u]"), 
       col = c("#4e6fc1", "#cd5cbf"), lty = c(1, 2), lwd = c(2, 1.5), 
       bty = "n", xpd = TRUE, cex = 0.9)

# Plot 3: Normalized population models comparison
plot(t_prime_post, N_prime_vals_post, type = "l", lwd = 2, col = "purple",
     xlab = "Normalized time (t')", ylab = "N' (normalized population)",
     main = "Normalized Population Dynamics: Model Comparison (post burn-in)",
     ylim = c(min(N_prime_vals_post) - 0.1, max(c(N_prime_vals_post, N_prime_null_vals_post, 
                      N_prime_slow_vals_post, K_prime_fast_vals_post)) * 1.1))
lines(t_prime_post, N_prime_null_vals_post, lwd = 2, col = "darkgray", lty = 3)
lines(t_prime_post, N_prime_slow_vals_post, lwd = 2, col = "#cd5cbf", lty = 3)
lines(t_prime_post, K_prime_fast_vals_post, lwd = 2, col = "#14eaa0", lty = 3)
# Adjusted legend position - higher and with less vertical height
legend(x = max(t_prime_post) * 1.02, y = max(c(N_prime_vals_post, N_prime_null_vals_post, 
                N_prime_slow_vals_post, K_prime_fast_vals_post)) * 1.05,
       legend = c("Dynamics", "Null", "Slow", "Fast"),
       col = c("purple", "darkgray", "#cd5cbf", "#14eaa0"),
       lty = c(1, 2, 3, 4), lwd = c(2, 2, 2, 2), bty = "n", xpd = TRUE, cex = 0.9,
       y.intersp = 0.8)  # Reduced spacing between legend items

# Add main title and subtitle with parameters and metrics
title(main = paste("Normalized Model with Burn-in Period"), 
      outer = TRUE, line = 4, cex.main = 1.2)
mtext(param_string, side = 3, outer = TRUE, line = 2.5, cex = 0.9)
mtext(metrics_string, side = 3, outer = TRUE, line = 1, cex = 0.8, col = "darkblue")