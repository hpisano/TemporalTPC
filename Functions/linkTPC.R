# ==================== SIMPLE FUNCTIONS =====================

# ==================== SIMPLE FUNCTIONS =====================

#' Simple growth rate function
#' 
#' @param T Temperature (either numeric constant or function T(t))
#' @param T_opt Optimal temperature in Celsius
#' @param epsilon Thermal sensitivity parameter (T_cmax - T_opt)
#' @param d_inf Rate at low temperature limit (d_{T-∞})
#' @param delta_r Scaling factor (Δ_r)
#' @return Either a numeric rate or a function r(t)
simple_growth_rate_r <- function(T, T_opt, epsilon, d_inf, delta_r) {
  # Check for valid epsilon
  if (epsilon <= 0) {
    stop("epsilon (T_cmax - T_opt) must be greater than 0")
  }
  
  # If T is a function, return a function of time
  if (is.function(T)) {
    return(function(t) {
      T_val <- T(t)
      x <- (T_val - T_opt) / epsilon
      rate_value <- exp(x) * (1 - x) - d_inf
      delta_r * rate_value
    })
  } else {
    # Otherwise, T is numeric, compute directly
    x <- (T - T_opt) / epsilon
    rate_value <- exp(x) * (1 - x) - d_inf
    return(delta_r * rate_value)
  }
}

#' u parametrized with the new 3 params
#' 
#' @param P_offset adimensional parameter for the offset between mu and Topt
#' @param P_amp adimensional amplitude
#' @param d_inf Rate at low temperature limit (d_{T-∞})
#' @param z Normalized temperature (either numeric constant or function z(t_prime))
#' @return Either a numeric rate or a function u(t_prime)
growth_rate_3dim_u <- function(P_offset, P_amp, d_inf, z) {
  # If z is a function, return a function of normalized time
  if (is.function(z)) {
    return(function(t_prime) {
      z_val <- z(t_prime)
      x <- P_offset + P_amp * z_val
      rate_value <- exp(x) * (1 - x) - d_inf
      return(rate_value)
    })
  } else {
    # Otherwise, z is numeric, return a constant function
    # that ignores its input (for compatibility with deSolve)
    x <- P_offset + P_amp * z
    rate_value <- exp(x) * (1 - x) - d_inf
    return(function(t_prime) rep(rate_value, length(t_prime)))
  }
}

# ==================== COMPLEX FUNCTIONS ====================

#' Complex Birth rate function based on UTPC
#' @param T Temperature in Celsius
#' @param T_opt Optimal temperature
#' @param epsilon Thermal sensitivity parameter
#' @return Birth rate
complex_birth_rate <- function(T, T_opt, epsilon) {
  x <- (T - T_opt) / epsilon
  return(exp(x) * (1 - x))
}

#' Simplified death rate function
#' @param T Temperature in Celsius
#' @param T_opt Optimal temperature
#' @param D Base death rate
#' @param epsilon Thermal sensitivity parameter
#' @param T0 Reference temperature in Kelvin
#' @param beta Cold mortality coefficient
#' @return Death rate
complex_death_rate <- function(T, T_opt, D, epsilon, T0, beta) {
  term1 <- (epsilon * T0) / (epsilon * T0 + beta)
  term2 <- exp((T - T_opt) / (epsilon * T0))
  term3 <- beta / (epsilon * T0 + beta)
  term4 <- exp((T_opt - T) / beta)
  
  return(D * (term1 * term2 + term3 * term4))
}

#' Growth rate function
#' @param T Temperature in Celsius
#' @param T_opt Optimal temperature
#' @param epsilon Thermal sensitivity parameter
#' @param D Base death rate
#' @param T0 Reference temperature in Kelvin
#' @param beta Cold mortality coefficient
#' @return Growth rate (birth - death)
complex_growth_rate <- function(T, T_opt, epsilon, D, T0, beta) {
  b <- complex_birth_rate(T, T_opt, epsilon)
  d <- complex_death_rate(T, T_opt, D, epsilon, T0, beta)
  return(b - d)
}

