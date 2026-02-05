#' Plot Birth, Death, and Growth Rates as Functions of Temperature
#' 
#' This script implements the Universal Thermal Performance Curve (UTPC) model
#' for birth rate and a simplified death rate model to calculate growth rate
#' as r(T) = b(T) - d(T).

# ==================== USER PARAMETERS ====================
# Adjust these values to modify the model behavior

# Temperature parameters (in degrees Celsius)
T_opt <- 25      # Optimal temperature for birth
TPCmax <- 28     # Max temp for birth
T0 <- 273     # Reference temperature in Kelvin (20°C = 293.15K)

# Death rate parameters
D <- 0.1         # Base death rate in permissive range
beta <- 6        # Cold tolerance coefficient
epsilon <- 0.008

# Permissive range boundaries (optional, for visualization)
T_cmin <- 15     # Lower limit of permissive range
T_cmax <- 29     # Upper limit of permissive range (not necessarily TPC_max)

# Plotting parameters
T_min <- 0       # Minimum temperature for plot
T_max <- 40      # Maximum temperature for plot
n_points <- 1000 # Number of points for smooth curves
ymin= -1.5
ymax= 1.5

# ==================== CALCULATIONS ====================

# Create temperature sequence
T_seq <- seq(T_min, T_max, length.out = n_points)

# Calculate rates
b_vals <- complex_birth_rate(T_seq, T_opt, epsilon)
d_vals <- complex_death_rate(T_seq, T_opt, D, epsilon, T0, beta)
r_vals <- b_vals - d_vals

# ==================== PLOTTING ====================

# Set up color scheme
colors <- c("Birth" = "blue", "Death" = "red", "Growth" = "green")

# Create main plot with fixed y-axis limits
par(mar = c(5, 4, 4, 8)) # Adjust margins for legend

plot(T_seq, b_vals, type = "l", col = colors["Birth"], lwd = 2,
     xlab = "Temperature (°C)", ylab = "Rate",
     main = "Thermal Rates Curve",
     ylim = c(ymin, ymax),  # Fixed y-axis limits
     xlim = c(T_min, T_max),
     cex.main = 1.2, cex.lab = 1.1, las = 1)

# Add other curves
lines(T_seq, d_vals, col = colors["Death"], lwd = 2, lty = 1)
lines(T_seq, r_vals, col = colors["Growth"], lwd = 2, lty = 1)

# Add horizontal line at y=0
abline(h = 0, lty = 3, col = "gray")

# Add horizontal reference lines at y = 10 and y = -10 (optional)
abline(h = 10, lty = 3, col = "lightgray", lwd = 0.5)
abline(h = -10, lty = 3, col = "lightgray", lwd = 0.5)

# Add vertical lines for critical temperatures
abline(v = T_opt, lty = 2, col = "darkgreen", lwd = 1.5)
abline(v = T_cmin, lty = 3, col = "blue", lwd = 1)
abline(v = T_cmax, lty = 3, col = "red", lwd = 1)

# Add legend
legend(par("usr")[2] + 0.5, par("usr")[4], 
       legend = c("Birth rate", 
                  "Death rate", 
                  "Growth rate",
                  paste0("T_opt = ", T_opt, "°C"),
                  paste0("T_cmin = ", T_cmin, "°C"),
                  paste0("T_cmax = ", T_cmax, "°C")),
       col = c(colors, "darkgreen", "blue", "red"),
       lty = c(1, 1, 1, 2, 3, 3),
       lwd = c(2, 2, 2, 1.5, 1, 1),
       cex = 0.8, bg = "white", xpd = TRUE)

# Add parameter information box
param_text <- paste(
  sprintf("Parameters:\nε = %.1f\nD = %.2f\nβ = %.1f\nT0 = %.1fK", 
          epsilon, D, beta, T0),
  sep = ""
)

# Place text in top left corner
text(T_min + 2, 9,  # Positioned at top left, just below y=10
     param_text, adj = c(0, 1), cex = 0.7, col = "darkgray")




