# =========================== POPULATION STATISTICS FUNCTIONS =================

#' Coefficient of Variation for population (continuous or discrete)
#' 
#' @param t_time Time vector (continuous, can be irregular)
#' @param N_pop Population vector (same length as t_time)
#' @param weighted Logical, if TRUE weights by time intervals for irregular sampling
#' @return Coefficient of variation CV = σ/μ
#' @details For irregular time series (common in ODE solutions), weights by time intervals
population_cv <- function(t_time, N_pop, weighted = TRUE) {
  if (length(t_time) != length(N_pop)) {
    stop("t_time and N_pop must have the same length")
  }
  
  if (any(N_pop < 0)) {
    warning("Negative population values detected. Results may be non-sensical.")
  }
  
  if (weighted && length(t_time) > 1) {
    # Weight by time intervals for irregular sampling
    dt <- diff(t_time)
    weights <- c(dt[1], dt)  # Extend first weight
    weights <- weights / sum(weights)
    
    # Weighted mean
    mean_N <- sum(N_pop * weights, na.rm = TRUE)
    
    # Weighted variance
    var_N <- sum(weights * (N_pop - mean_N)^2, na.rm = TRUE)
    sd_N <- sqrt(var_N)
  } else {
    # Regular sampling (equal weights)
    mean_N <- mean(N_pop, na.rm = TRUE)
    sd_N <- sd(N_pop, na.rm = TRUE)
  }
  
  if (abs(mean_N) < .Machine$double.eps) {
    return(Inf)  # Division by zero
  }
  
  return(sd_N / mean_N)
}

#' Autocorrelation Function for population time series
#' 
#' @param t_time Time vector
#' @param N_pop Population vector (same length as t_time)
#' @param tau_max Maximum lag to compute (defaults to 1/4 of time range)
#' @param n_lags Number of lag points to compute (default 50)
#' @param method "fft" for fast computation on regular grids, "direct" for irregular
#' @return List with lags and corresponding autocorrelation values
#' @details Efficient ACF computation suitable for ODE outputs
population_acf <- function(t_time, N_pop, tau_max = NULL, n_lags = 50, method = "fft") {
  if (length(t_time) != length(N_pop)) {
    stop("t_time and N_pop must have the same length")
  }
  
  # Remove any NA values
  valid <- !is.na(t_time) & !is.na(N_pop)
  t_time <- t_time[valid]
  N_pop <- N_pop[valid]
  
  if (length(t_time) < 2) {
    stop("Need at least 2 valid time points")
  }
  
  # Detrend by subtracting mean
  N_detrend <- N_pop - mean(N_pop, na.rm = TRUE)
  
  # Check if time is regularly sampled
  dt <- diff(t_time)
  is_regular <- all(abs(dt - mean(dt)) < 1e-10 * mean(dt))
  
  if (is_regular && method == "fft") {
    # Use FFT for regular sampling (much faster)
    return(acf_regular(t_time, N_detrend, tau_max, n_lags))
  } else {
    # Use interpolation for irregular or direct computation
    return(acf_irregular(t_time, N_detrend, tau_max, n_lags))
  }
}

#' Helper: ACF for regularly sampled data using FFT
acf_regular <- function(t_time, N_detrend, tau_max = NULL, n_lags = 50) {
  dt <- t_time[2] - t_time[1]
  n <- length(N_detrend)
  
  if (is.null(tau_max)) {
    tau_max <- (t_time[n] - t_time[1]) / 4
  }
  
  # Compute ACF using FFT (much faster for large datasets)
  N_fft <- fft(N_detrend)
  acf_raw <- Re(fft(Conj(N_fft) * N_fft, inverse = TRUE)) / (n * var(N_detrend))
  acf_raw <- acf_raw[1:n] / acf_raw[1]  # Normalize
  
  # Create lags
  max_lag_idx <- min(n, round(tau_max / dt))
  lags <- seq(0, tau_max, length.out = n_lags)
  lag_idxs <- round(lags / dt) + 1
  lag_idxs <- pmin(lag_idxs, n)  # Don't exceed array bounds
  
  acf_vals <- acf_raw[lag_idxs]
  
  return(list(lags = lags, acf = acf_vals, method = "fft"))
}

#' Helper: ACF for irregularly sampled data using interpolation
acf_irregular <- function(t_time, N_detrend, tau_max = NULL, n_lags = 50) {
  t_range <- range(t_time)
  t_total <- t_range[2] - t_range[1]
  
  if (is.null(tau_max)) {
    tau_max <- t_total / 4
  }
  
  # Create interpolation function
  N_interp <- approxfun(t_time, N_detrend, rule = 2)
  
  # Create lags
  lags <- seq(0, tau_max, length.out = n_lags)
  
  # Function to compute ACF for a given lag
  compute_acf_for_lag <- function(tau) {
    # Only consider times where t + tau is within range
    valid_t <- t_time[t_time <= t_range[2] - tau]
    
    if (length(valid_t) == 0) {
      return(0)
    }
    
    # Compute product at valid times
    products <- N_interp(valid_t) * N_interp(valid_t + tau)
    
    # Average the products
    mean_product <- mean(products, na.rm = TRUE)
    
    # Normalize by variance at lag 0
    variance <- var(N_detrend, na.rm = TRUE)
    
    if (variance < .Machine$double.eps) {
      return(0)
    }
    
    return(mean_product / variance)
  }
  
  # Compute ACF for all lags
  acf_vals <- sapply(lags, compute_acf_for_lag)
  
  return(list(lags = lags, acf = acf_vals, method = "interpolation"))
}

# =========================== COMPARISON FUNCTIONS ===========================

#' Calculate absolute deviation between two functions
#' @param t_f Final time for integration
#' @param N_func Function N(t) (can be numeric vector or function)
#' @param m_func Function m(t) (can be numeric vector or function)
#' @return Absolute deviation Dev_m(t_f) = ∫₀^{t_f} |N(t) - m(t)| dt
#' @details Computes the absolute deviation between two functions over [0, t_f]
absolute_deviation <- function(N_func, m_func, t_f, tol = 1e-4) {
  integrand <- function(t) {
    abs(N_func(t) - m_func(t))
  }
  
  # Handle vector input
  integrand_vector <- function(x) {
    sapply(x, integrand)
  }
  
  result <- cubature::adaptIntegrate(integrand_vector, lowerLimit = 0, upperLimit = t_f,
                           tol = tol, fDim = 1)
  result$integral
}

