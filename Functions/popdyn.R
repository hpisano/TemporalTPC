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

create_normalized_systems <- function(u_func, P_time, N_prime0, threshold = 1e-2) {
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
  
  # Helper function to apply dynamic threshold in ODE
  dNprime_dynamic <- function(N_prime, u_val) {
    # Apply threshold: if N_prime is below threshold, set derivative to 0
    if (N_prime < threshold) {
      return(0)
    }
    return(P_time * (u_val * N_prime - N_prime^2))
  }
  
  list(
    dynamic = function(t, y, params) {
      N_prime <- y[1]
      u_val <- u_func(t)
      dNprime <- dNprime_dynamic(N_prime, u_val)
      list(dNprime)
    },
    
    null = function(t, y, params) {
      N_prime <- y[1]
      u_val <- u_null_func(t)
      dNprime <- dNprime_dynamic(N_prime, u_val)
      list(dNprime)
    },
    
    slow = function(t, y, params) {
      N_prime <- y[1]
      u_val <- u_slow_func(t)
      dNprime <- dNprime_dynamic(N_prime, u_val)
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
  
  # Calculate the integral between t=0 and t=1 of P_time * u_func(t)
  integral_result <- tryCatch({
    integrand_func <- function(t) {
      return(P_time * u_func(t))
    }
    
    result <- integrate(integrand_func, lower = 0, upper = 1)$value
    
    if (!is.finite(result)) {
      t_vals <- seq(0, 1, length.out = 1000)
      integrand_vals <- sapply(t_vals, integrand_func)
      result <- mean(integrand_vals) * (1 - 0)
    }
    
    result
  }, error = function(e) {
    NA_real_
  })
  
  # Helper function to solve ODE with multiple fallback methods
  solve_ode_robust <- function(system_func, y0, times, max_attempts = 4) {
    methods_to_try <- c("lsoda", "bdf", "radau", "rk4")
    
    for (method in methods_to_try[1:max_attempts]) {
      solution <- tryCatch({
        # Suppress the DLSODA warnings
        suppressWarnings({
          if (method %in% c("bdf", "radau")) {
            # These methods are better for stiff systems
            ode(y = y0, 
                times = times, 
                func = system_func, 
                parms = NULL,
                method = method,
                atol = 1e-6, rtol = 1e-6,
                maxsteps = 50000)
          } else if (method == "rk4") {
            # Fixed step method as last resort
            ode(y = y0, 
                times = times, 
                func = system_func, 
                parms = NULL,
                method = "rk4")
          } else {
            ode(y = y0, 
                times = times, 
                func = system_func, 
                parms = NULL,
                method = method,
                atol = 1e-6, rtol = 1e-6,
                maxsteps = 50000)
          }
        })
      }, error = function(e) NULL, warning = function(w) NULL)
      
      # Check if solution is valid
      if (!is.null(solution) && 
          nrow(solution) == length(times) && 
          all(is.finite(solution[, "N_prime"]))) {
        return(solution)
      }
    }
    
    # If all methods fail, return NULL
    return(NULL)
  }
  
  # Try to solve ODEs with error handling
  tryCatch({
    # Solve ODE for dynamic system
    solution <- solve_ode_robust(pop_systems_prime$dynamic, 
                                  c(N_prime = N_prime0), 
                                  t_prime)
    
    # Solve null system
    solution_null <- solve_ode_robust(pop_systems_prime$null, 
                                       c(N_prime = N_prime0), 
                                       t_prime)
    
    # Solve slow system
    solution_slow <- solve_ode_robust(pop_systems_prime$slow, 
                                       c(N_prime = N_prime0), 
                                       t_prime)
    
    # Check if any ODE solutions failed - return NA instead of stopping
    if (is.null(solution) || is.null(solution_null) || is.null(solution_slow)) {
      return(data.frame(
        P_offset = P_offset,
        P_amp = P_amp,
        P_time = P_time,
        E_N_prime = NA_real_,
        speed = NA_real_,
        abs_dev_slow = NA_real_,
        abs_dev_fast = NA_real_,
        integral_P_time_u = integral_result
      ))
    }
    
    # Apply threshold to prevent extremely small population values
    threshold_clean <- function(mat, thresh = 0.01) {
      if (is.null(mat) || nrow(mat) == 0) return(mat)
      
      nprime_col <- which(colnames(mat) == "N_prime")
      if (length(nprime_col) == 0) nprime_col <- 2
      
      finite_mask <- is.finite(mat[, nprime_col])
      below_thresh <- mat[, nprime_col] < thresh & finite_mask
      mat[below_thresh, nprime_col] <- 0
      
      return(mat)
    }
    
    solution <- threshold_clean(solution)
    solution_null <- threshold_clean(solution_null)
    solution_slow <- threshold_clean(solution_slow)
    
    # Extract post-burn-in values
    post_burn_idx <- solution[, "time"] > burn_in_time
    N_prime_post_burn <- solution[post_burn_idx, "N_prime"]
    
    # Check if we have valid post-burn-in data - return NA instead of stopping
    if (length(N_prime_post_burn) == 0 || all(!is.finite(N_prime_post_burn))) {
      return(data.frame(
        P_offset = P_offset,
        P_amp = P_amp,
        P_time = P_time,
        E_N_prime = NA_real_,
        speed = NA_real_,
        abs_dev_slow = NA_real_,
        abs_dev_fast = NA_real_,
        integral_P_time_u = integral_result
      ))
    }
    
    # Calculate expected N_prime after burn-in (ignore NAs)
    E_N_prime <- mean(N_prime_post_burn[is.finite(N_prime_post_burn)], na.rm = TRUE)
    
    if (!is.finite(E_N_prime)) {
      E_N_prime <- NA_real_
    }
    
    # Create functions for each system
    N_prime_dynamic_func <- approxfun(solution[, "time"], solution[, "N_prime"], rule = 2)
    N_prime_null_func <- approxfun(solution_null[, "time"], solution_null[, "N_prime"], rule = 2)
    N_prime_slow_func <- approxfun(solution_slow[, "time"], solution_slow[, "N_prime"], rule = 2)
    
    K_prime_fast_func <- pop_systems_prime$fast
    
    # Calculate deviations with error handling
    abs_dev_slow_prime <- tryCatch({
      absolute_deviation(t_f = max(t_prime),
                         t_0 = burn_in_time,
                         N_func = N_prime_dynamic_func,
                         m_func = N_prime_slow_func)
    }, error = function(e) NA_real_)
    
    abs_dev_fast_prime <- tryCatch({
      absolute_deviation(t_f = max(t_prime),
                         t_0 = burn_in_time,
                         N_func = N_prime_dynamic_func,
                         m_func = K_prime_fast_func)
    }, error = function(e) NA_real_)
    
    if (!is.finite(abs_dev_slow_prime)) abs_dev_slow_prime <- NA_real_
    if (!is.finite(abs_dev_fast_prime)) abs_dev_fast_prime <- NA_real_
    
    # Calculate speed metric
    speed <- NA_real_
    if (is.finite(abs_dev_fast_prime) && is.finite(abs_dev_slow_prime)) {
      denominator <- abs_dev_fast_prime + abs_dev_slow_prime
      if (denominator > 0) {
        speed <- abs_dev_slow_prime / denominator
      } else {
        speed <- 0
      }
    }
    
    if (is.finite(speed)) {
      speed <- max(0, min(1, speed))
    }
    
    # Return results
    return(data.frame(
      P_offset = P_offset,
      P_amp = P_amp,
      P_time = P_time,
      E_N_prime = E_N_prime,
      speed = speed,
      abs_dev_slow = abs_dev_slow_prime,
      abs_dev_fast = abs_dev_fast_prime,
      integral_P_time_u = integral_result
    ))
    
  }, error = function(e) {
    # Return NA values for failed simulations (silently)
    return(data.frame(
      P_offset = P_offset,
      P_amp = P_amp,
      P_time = P_time,
      E_N_prime = NA_real_,
      speed = NA_real_,
      abs_dev_slow = NA_real_,
      abs_dev_fast = NA_real_,
      integral_P_time_u = integral_result
    ))
  })
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

