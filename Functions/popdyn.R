# ==================== POPULATION GROWTH FUNCTIONS =====================

#' Logistic growth model
#' 
#' @param N Population size (numeric or function N(t))
#' @param r Intrinsic growth rate (numeric or function r(t))
#' @param alpha Density-dependent mortality coefficient (numeric or function alpha(t))
#' @return Rate of change dN/dt (numeric or function dN/dt(t))
logistic_growth <- function(N, r, alpha) {
  # Determine if we need to return a function
  is_dynamic <- is.function(N) || is.function(r) || is.function(alpha)
  
  if (is_dynamic) {
    # Return a function that can be evaluated at any time t
    return(function(t) {
      # Evaluate N at time t (if it's a function)
      N_val <- if (is.function(N)) N(t) else N
      
      # Evaluate r at time t (if it's a function)
      r_val <- if (is.function(r)) r(t) else r
      
      # Evaluate alpha at time t (if it's a function)
      alpha_val <- if (is.function(alpha)) alpha(t) else alpha
      
      # Compute and return the rate
      return((r_val - alpha_val * N_val) * N_val)
    })
  } else {
    # All inputs are numeric, compute directly
    return((r - alpha * N) * N)
  }
}

#' 3-parameter adimensional population growth model
#' 
#' @param N_prime Normalized population size (numeric or function N'(t'))
#' @param P_time Scaled time parameter (numeric or function P_time(t'))
#' @param u Thermal response function (numeric or function u(t'))
#' @return Rate of change dN'/dt' (numeric or function dN'/dt'(t'))
normalized_3dim_growth <- function(N_prime, P_time, u) {
  
  # Check if we need dynamic evaluation
  is_dynamic <- is.function(N_prime) || is.function(P_time) || is.function(u)
  
  if (is_dynamic) {
    # Create a function for dynamic evaluation
    return(function(t_prime) {
      # Handle vector input
      if (length(t_prime) > 1) {
        return(sapply(t_prime, function(tp) {
          N_val <- if (is.function(N_prime)) N_prime(tp) else N_prime
          P_val <- if (is.function(P_time)) P_time(tp) else P_time
          u_val <- if (is.function(u)) u(tp) else u
          return(P_val * (u_val * N_val - N_val^2))
        }))
      } else {
        # Scalar input
        N_val <- if (is.function(N_prime)) N_prime(t_prime) else N_prime
        P_val <- if (is.function(P_time)) P_time(t_prime) else P_time
        u_val <- if (is.function(u)) u(t_prime) else u
        return(P_val * (u_val * N_val - N_val^2))
      }
    })
  } else {
    # Static computation
    return(P_time * (u * N_prime - N_prime^2))
  }
}

#' Create continuous ODE system for all population models
#' 
#' @param r_func Continuous function r(t) for dynamic model
#' @param alpha Density-dependent mortality coefficient
#' @param N0 Initial population size
#' @return List of ODE functions for dynamic, null, slow, and fast models
create_population_systems <- function(r_func, alpha, N0) {
  # Extract mean r value once (for slow model)
  # We need to evaluate over a reasonable time range to get mean
  # Let's use a fixed time range for mean calculation
  t_range <- c(0, 100)  # Adjust as needed
  t_eval <- seq(t_range[1], t_range[2], length.out = 1000)
  r_vals <- r_func(t_eval)
  r_mean <- mean(r_vals)
  
  # Create constant functions for null and slow models
  r_null_func <- function(t) r_func(0)  # Initial r value (or any constant)
  r_slow_func <- function(t) r_mean
  
  # Return ODE functions for each model
  list(
    # Dynamic model ODE
    dynamic = function(t, y, params) {
      r <- r_func(t)
      N <- y[1]
      dN <- (r - alpha * N) * N
      list(dN)
    },
    
    # Null model ODE (constant r at initial value)
    null = function(t, y, params) {
      r <- r_null_func(t)
      N <- y[1]
      dN <- (r - alpha * N) * N
      list(dN)
    },
    
    # Slow model ODE (constant r at mean value)
    slow = function(t, y, params) {
      r <- r_slow_func(t)
      N <- y[1]
      dN <- (r - alpha * N) * N
      list(dN)
    },
    
    # Fast model is just r(t)/alpha (carrying capacity)
    fast = function(t) {
      r_func(t) / alpha
    }
  )
}

create_normalized_systems <- function(u_func, P_time, N_prime0) {
  # Calculate mean u for slow model
  t_range <- c(0, 10)
  t_eval <- seq(t_range[1], t_range[2], length.out = 1000)
  u_vals <- u_func(t_eval)
  u_mean <- mean(u_vals)
  
  # Get initial value of u for null model
  u_initial <- u_func(0)
  
  # Constant functions for null and slow models
  u_null_func <- function(t) u_initial  # Use initial u value, not 0
  u_slow_func <- function(t) u_mean
  
  list(
    dynamic = function(t, y, params) {
      N_prime <- y[1]
      u_val <- u_func(t)
      dNprime <- P_time * (u_val * N_prime - N_prime^2)
      list(dNprime)
    },
    
    null = function(t, y, params) {
      N_prime <- y[1]
      u_val <- u_null_func(t)
      dNprime <- P_time * (u_val * N_prime - N_prime^2)
      list(dNprime)
    },
    
    slow = function(t, y, params) {
      N_prime <- y[1]
      u_val <- u_slow_func(t)
      dNprime <- P_time * (u_val * N_prime - N_prime^2)
      list(dNprime)
    },
    
    fast = function(t) {
      u_func(t)
    }
  )
}


run_range_3dim_simulation <- function(P_offset, P_amp, P_time) {
  # Use the global parameters defined at the beginning
  z_func <- create_z_function(model = z_params$model,
                              period = z_params$period,
                              phase = z_params$phase)
  
  u_func <- growth_rate_3dim_u(P_offset = P_offset,
                               P_amp = P_amp,
                               d_inf = d_inf,
                               z = z_func)
  
  pop_systems_prime <- create_normalized_systems(u_func, P_time, N_prime0)
  
  # Solve ODE for dynamic system
  solution <- ode(y = c(N_prime = N_prime0), 
                  times = t_prime, 
                  func = pop_systems_prime$dynamic, 
                  parms = NULL, method = "vode", mf = 22)
  
  # Solve null system
  solution_null <- ode(y = c(N_prime = N_prime0), 
                       times = t_prime, 
                       func = pop_systems_prime$null, 
                       parms = NULL)
  
  # Solve slow system
  solution_slow <- ode(y = c(N_prime = N_prime0), 
                       times = t_prime, 
                       func = pop_systems_prime$slow, 
                       parms = NULL)
  
  # Apply threshold to prevent extremely small population values
  solution[, "N_prime"][solution[, "N_prime"] < 0.01] <- 0
  solution_null[, "N_prime"][solution_null[, "N_prime"] < 0.01] <- 0
  solution_slow[, "N_prime"][solution_slow[, "N_prime"] < 0.01] <- 0
  
  # Extract post-burn-in values
  post_burn_idx <- solution[, "time"] > burn_in_time
  N_prime_post_burn <- solution[post_burn_idx, "N_prime"]
  
  # Create functions for each system
  N_prime_dynamic_func <- approxfun(solution[, "time"], solution[, "N_prime"], rule = 2)
  N_prime_null_func <- approxfun(solution_null[, "time"], solution_null[, "N_prime"], rule = 2)
  N_prime_slow_func <- approxfun(solution_slow[, "time"], solution_slow[, "N_prime"], rule = 2)
  
  # Fast system function is u_func(t) as defined in create_normalized_systems
  K_prime_fast_func <- pop_systems_prime$fast
  
  # Calculate deviations
  abs_dev_slow_prime <- absolute_deviation(t_f = max(t_prime),
                                           t_0 = burn_in_time,
                                           N_func = N_prime_dynamic_func,
                                           m_func = N_prime_slow_func)
  
  abs_dev_fast_prime <- absolute_deviation(t_f = max(t_prime),
                                           t_0 = burn_in_time,
                                           N_func = N_prime_dynamic_func,
                                           m_func = K_prime_fast_func)
  
  # Calculate speed metric
  if (abs_dev_fast_prime + abs_dev_slow_prime == 0) {
    speed <- 0
  } else {
    speed <- abs_dev_slow_prime / (abs_dev_fast_prime + abs_dev_slow_prime)
  }
  
  # Return all metrics including the absolute deviations
  return(data.frame(
    P_offset = P_offset,
    P_amp = P_amp,
    P_time = P_time,
    E_N_prime = mean(N_prime_post_burn, na.rm = TRUE),
    speed = speed,
    abs_dev_slow = abs_dev_slow_prime,
    abs_dev_fast = abs_dev_fast_prime,
    stringsAsFactors = FALSE
  ))
}

# Function to create heatmap for a specific P_offset value
create_heatmap <- function(data, metric, title_prefix, color_palette = "viridis") {
  # Determine appropriate limits based on metric
  if (metric == "E_N_prime") {
    fill_limits <- c(0, 1)
    fill_label <- "E[N']"
  } else if (metric == "speed") {
    fill_limits <- c(0, 1)
    fill_label <- "Speed"
  } else {
    fill_limits <- NULL
    fill_label <- metric
  }
  
  ggplot(data, aes(x = P_amp, y = P_time, fill = !!sym(metric))) +
    geom_tile() +
    scale_fill_viridis(option = color_palette, limits = fill_limits) +
    labs(
      title = paste(title_prefix, "\nP_offset =", unique(data$P_offset)),
      x = "P_amp",
      y = "P_time",
      fill = fill_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      legend.position = "right"
    ) +
    coord_fixed(ratio = (max(data$P_amp) - min(data$P_amp)) / 
                  (max(data$P_time) - min(data$P_time)) * 0.8)
}

