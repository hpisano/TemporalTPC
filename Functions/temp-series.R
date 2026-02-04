# ==================== T-FUNCTIONS =====================

#' Generate deterministic temperature component as a function
#' 
#' @param model Character specifying the model type
#' @param ... Parameters for the selected model
#' @return A function T(t) that can be evaluated at any t
generate_deterministic_T <- function(model = "linear", ...) {
  params <- list(...)
  
  switch(model,
    "linear" = {
      if (is.null(params$a) || is.null(params$b)) {
        stop("For 'linear' model, provide a and b")
      }
      a <- params$a
      b <- params$b
      function(t) a * t + b
    },
    
    "sine" = {
      if (is.null(params$mu) || is.null(params$sigma) || is.null(params$p)) {
        stop("For 'sine' model, provide mu, sigma, and p")
      }
      mu <- params$mu
      sigma <- params$sigma
      p <- params$p
      function(t) mu + sigma * sin(2 * pi * t / p)
    },
    
    "sum_sines" = {
      if (is.null(params$mu_vec) || is.null(params$sigma_vec) || is.null(params$p_vec)) {
        stop("For 'sum_sines' model, provide mu_vec, sigma_vec, and p_vec")
      }
      if (!all(length(params$mu_vec) == length(params$sigma_vec), 
               length(params$sigma_vec) == length(params$p_vec))) {
        stop("mu_vec, sigma_vec, and p_vec must have the same length")
      }
      mu_vec <- params$mu_vec
      sigma_vec <- params$sigma_vec
      p_vec <- params$p_vec
      function(t) {
        result <- rep(0, length(t))
        for (i in seq_along(mu_vec)) {
          result <- result + mu_vec[i] + sigma_vec[i] * sin(2 * pi * t / p_vec[i])
        }
        result
      }
    },
    
    "sine_lin_amplitude" = {
      if (is.null(params$mu) || is.null(params$sigma) || is.null(params$p)) {
        stop("For 'sine_lin_amplitude' model, provide mu, sigma, and p")
      }
      mu <- params$mu
      sigma <- params$sigma
      p <- params$p
      function(t) mu + sigma * t * sin(2 * pi * t / p)
    },
    
    "sine_lin_mean" = {
      if (is.null(params$mu) || is.null(params$sigma) || is.null(params$p)) {
        stop("For 'sine_lin_mean' model, provide mu, sigma, and p")
      }
      mu <- params$mu
      sigma <- params$sigma
      p <- params$p
      function(t) mu * t + sigma * sin(2 * pi * t / p)
    },
    
    stop(paste("Unknown model type:", model))
  )
}

#' Generate Ornstein-Uhlenbeck process function
#' 
#' @param mu Long-term mean
#' @param T0 Initial value
#' @param theta Mean reversion rate
#' @param sigma Volatility
#' @param seed Random seed (optional)
#' @return A function OU(t) that generates the process at any t
generate_ou_process <- function(mu, T0, theta, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Store the last value and time for incremental updates
  last_t <- 0
  last_x <- T0
  
  function(t) {
    # Handle vector input
    if (length(t) > 1) {
      # Sort times and generate incrementally
      t_sorted <- sort(t)
      x <- numeric(length(t_sorted))
      x[1] <- last_x
      
      for (i in 2:length(t_sorted)) {
        dt <- t_sorted[i] - t_sorted[i-1]
        dW <- rnorm(1, mean = 0, sd = sqrt(dt))
        x[i] <- x[i-1] + theta * (mu - x[i-1]) * dt + sigma * dW
      }
      
      # Store the last values
      last_t <<- t_sorted[length(t_sorted)]
      last_x <<- x[length(x)]
      
      # Return in original order
      return(x[order(order(t))])
    } else {
      # Single time point
      dt <- t - last_t
      dW <- rnorm(1, mean = 0, sd = sqrt(abs(dt)))
      last_x <<- last_x + theta * (mu - last_x) * dt + sigma * dW
      last_t <<- t
      return(last_x)
    }
  }
}

#' Generate stochastic noise function
#' 
#' @param noise_type Type of noise: "ou", "white", or "none"
#' @param std Standard deviation for white noise
#' @param theta Mean reversion rate for OU noise
#' @param mu Mean for OU noise
#' @param seed Random seed (optional)
#' @return A function noise(t) that generates noise at any t
generate_noise <- function(noise_type = "none", std = 1, 
                          theta = 0.1, mu = 0, seed = NULL) {
  if (noise_type == "none") {
    return(function(t) rep(0, length(t)))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  switch(noise_type,
    "white" = {
      # For white noise, we need to generate new random values each call
      # but consistent for the same t values
      white_noise_cache <- new.env(hash = TRUE)
      function(t) {
        key <- paste(sort(t), collapse = ",")
        if (!exists(key, envir = white_noise_cache)) {
          assign(key, rnorm(length(t), mean = 0, sd = std), envir = white_noise_cache)
        }
        get(key, envir = white_noise_cache)[order(order(t))]
      }
    },
    
    "ou" = {
      # Create an OU process function
      ou_func <- generate_ou_process(mu = mu, T0 = 0, 
                                     theta = theta, sigma = std, 
                                     seed = seed)
      function(t) {
        noise <- ou_func(t)
        # Center the noise (optional)
        noise - mean(noise)
      }
    },
    
    stop("noise_type must be 'white', 'ou', or 'none'")
  )
}

#' Main temperature generation function
#' 
#' @param model Deterministic model type
#' @param noise_type Stochastic noise type
#' @param det_params List of deterministic model parameters
#' @param noise_params List of noise parameters
#' @param seed Random seed for reproducibility
#' @return A function T(t) that can be evaluated at any t
generate_T <- function(model = "linear",
                      noise_type = "none",
                      det_params = list(),
                      noise_params = list(std = 1, theta = 0.1, mu = 0),
                      seed = NULL) {
  
  # Generate deterministic component function
  det_func <- do.call(generate_deterministic_T, 
                     c(list(model = model), det_params))
  
  # Generate noise function if requested
  if (noise_type != "none") {
    noise_func <- do.call(generate_noise,
                         c(list(noise_type = noise_type, seed = seed), 
                           noise_params))
    function(t) det_func(t) + noise_func(t)
  } else {
    det_func
  }
}

#' Main temperature generation function
#' 
#' @param model Deterministic model type
#' @param noise_type Stochastic noise type
#' @param det_params List of deterministic model parameters
#' @param noise_params List of noise parameters
#' @param seed Random seed for reproducibility
#' @return A function T(t) that can be evaluated at any t
generate_T <- function(model = "linear",
                      noise_type = "none",
                      det_params = list(),
                      noise_params = list(std = 1, theta = 0.1, mu = 0),
                      seed = NULL) {
  
  # Generate deterministic component function
  det_func <- do.call(generate_deterministic_T, 
                     c(list(model = model), det_params))
  
  # Generate noise function if requested
  if (noise_type != "none") {
    noise_func <- do.call(generate_noise,
                         c(list(noise_type = noise_type, seed = seed), 
                           noise_params))
    function(t) det_func(t) + noise_func(t)
  } else {
    det_func
  }
}

# ====================  Z-FUNCTIONS =====================

#' Generate standardized function (mean=0, variance=1)
#' 
#' @param model Type of function: "sine", "ar1", "custom", or "from_temp"
#' @param ... Parameters for the selected model
#' @return A function z(t) with mean ≈ 0 and variance ≈ 1
create_z_function <- function(model = "sine", ...) {
  params <- list(...)
  
  switch(model,
    "sine" = {
      phase <- params$phase %||% 0
      amplitude <- sqrt(2)  # Ensures variance = 1
      function(t) amplitude * sin(2 * pi * t + phase)
    },
    
    "ar1" = {
      # For AR(1), we need to generate a continuous function
      # We can use an OU process as a continuous analog
      phi <- params$phi %||% 0.5
      seed <- params$seed
      
      # Create an OU process with parameters that match AR(1) properties
      # For continuous-time AR(1), the relationship is: phi = exp(-theta * delta_t)
      # We'll use delta_t = 1 for unit time steps
      theta <- -log(phi)  # Mean reversion rate
      sigma <- sqrt(1 - phi^2)  # Volatility to achieve unit variance
      
      ou_func <- generate_ou_process(mu = 0, T0 = 0, 
                                     theta = theta, sigma = sigma, 
                                     seed = seed)
      
      # Standardize the output
      function(t) {
        # We'll generate and standardize on the fly for the given t
        x <- ou_func(t)
        (x - mean(x)) / sd(x)
      }
    },
    
    "custom" = {
      t_values <- params$t_values
      values <- params$values
      
      if (is.null(t_values) || is.null(values)) {
        stop("For 'custom' model, provide t_values and values")
      }
      
      # Standardize
      values_std <- (values - mean(values)) / sd(values)
      
      # Create a proper interpolation function
      approxfun(t_values, values_std, method = "linear", rule = 2)
    },
    
    "from_temp" = {
      # This one needs a temperature function as input
      temp_func <- params$temp_func
      t_range <- params$t_range %||% c(0, 100)
      
      if (is.null(temp_func)) {
        stop("For 'from_temp' model, provide temp_func")
      }
      
      # Evaluate on a grid to get mean and sd
      t_grid <- seq(t_range[1], t_range[2], length.out = 1000)
      temp_vals <- temp_func(t_grid)
      mu <- mean(temp_vals)
      sigma <- sd(temp_vals)
      
      # Create standardized function
      function(t) (temp_func(t) - mu) / sigma
    },
    
    stop("Unknown model type. Use 'sine', 'ar1', 'custom', or 'from_temp'")
  )
}

#' Helper function: NULL coalescing operator
#' @param x First value
#' @param y Default value
#' @return x if not NULL, else y
`%||%` <- function(x, y) if (!is.null(x)) x else y