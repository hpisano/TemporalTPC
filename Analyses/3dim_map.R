# ==================== PARAMETERS =====================
library(deSolve, cubature)

# Normalized time settings
t_prime_start <- 0
t_prime_end <- 20
t_prime_step <- 0.01
t_prime <- seq(t_prime_start, t_prime_end, by = t_prime_step)

# ADD THIS: Burn-in period (in normalized time units)
burn_in_time <- 10  # Adjust this value as needed

# z (normalized temperature) parameters
z_params <- list(
  model = "sine",
  phase = 0
)

# u (thermal response) parameters
u_params <- list(
  P_offset = -4,
  P_amp = 0.5,
  d_inf = 0.2
)

# Normalized population parameters
N_prime_params <- list(
  N_prime0 = 0.5,
  P_time = 3
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
# (or use the original functions but integrate from burn_in_time)
abs_dev_null_prime <- absolute_deviation(t_f = max(t_prime),
                                        t_0 = burn_in_time,  # ADD THIS if your function supports it
                                        N_func = N_prime_dynamic_func,
                                        m_func = N_prime_null_func)

abs_dev_slow_prime <- absolute_deviation(t_f = max(t_prime),
                                        t_0 = burn_in_time,  # ADD THIS
                                        N_func = N_prime_dynamic_func,
                                        m_func = N_prime_slow_func)

abs_dev_fast_prime <- absolute_deviation(t_f = max(t_prime),
                                        t_0 = burn_in_time,  # ADD THIS
                                        N_func = N_prime_dynamic_func,
                                        m_func = K_prime_fast_func)

# If your absolute_deviation function doesn't support t_0 parameter, 
# you may need to modify it or create a wrapper:
# absolute_deviation_burn <- function(t_0, t_f, N_func, m_func, ...) {
#   integrand <- function(tau) {
#     abs(N_func(tau) - m_func(tau))
#   }
#   integrate(integrand, lower = t_0, upper = t_f, ...)$value
# }

speed = abs_dev_slow_prime / (abs_dev_fast_prime + abs_dev_slow_prime)

# ==================== PLOTTING (WITH BURN-IN) =====================

