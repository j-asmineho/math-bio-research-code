
##### Allow to settle after times 1000, then plot only for the last 50 years
library(deSolve)

# Parameters
N <- 5000000                      # Total population
mu <- 1 / 50                      # Birth/death rate (annual)
gamma <- 365 / 13                 # Recovery rate (13-day infectious period)
R0 <- 17                          # Basic reproduction number
beta <- R0 * (gamma + mu)         # Transmission rate
alpha <- 0.1                      # Seasonal forcing amplitude
p <- NA                           # Placeholder for vaccination rate

# Combine into a parameter vector
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p
)

# Equilibrium
Seqm <- (with(as.list(parms), N / R0))
Ieqm <- (with(as.list(parms), N * mu / (gamma + mu) * (1 - 1 / R0)))

# Initial state near equilibrium
Si <- round((Seqm + (N-Seqm)*0.01))  # Slightly perturbed from equilibrium
Ii <- round((Ieqm * 0.9))                # 90% of equilibrium infectious population
Ri <- round(N - Si - Ii)                 # Remaining population in recovered state
ic <- c(S = Si, I = Ii, R = Ri)          # Initial condition vector

# SIR Model with Seasonal Forcing and Vaccination
SIR.vector.field <- function(t, vars, parms) {
  with(as.list(c(parms, vars)), {
    # Seasonal forcing for transmission rate
    beta_t <- beta * (1 + alpha * cos(2 * pi * t))
    
    # Differential equations
    dS <- mu * (1 - p) * N - beta_t * S * I / N - mu * S
    dI <- beta_t * S * I / N - gamma * I - mu * I
    dR <- gamma * I + mu * p * N - mu * R
    
    # Return the derivatives
    list(c(dS, dI, dR))
  })
}

# Periodogram Function
periodogram_d <- function(df, xlim = c(0, 10), color, ...) {
  with(df, {
    s <- spectrum(I, plot = FALSE)
    adjusted_freq <- 52 * s$freq  # Convert from weeks to years
    periods <- 1 / adjusted_freq
    plot(
      periods, s$spec, type = "l", col = color, lwd = 2,
      xlab = "Period (Years)", ylab = "Spectral Density",
      xlim = xlim, ...
    )
  })
}

# Define a sequence of p values to iterate over
p_values <- c(0, 0.3, 0.45, 0.66, 0.7)
colors <- rainbow(length(p_values))

# Time sequence for the simulation
times <- seq(0,500, by = 1/52)

par(mfrow = c(length(p_values), 2), mar = c(4, 4, 2, 1))

# Loop over each value of p
for (i in seq_along(p_values)) {
  p <- p_values[i]
  
  # Update parameters
  parms <- c(
    beta = beta,
    gamma = gamma,
    mu = mu,
    alpha = alpha,
    N = N,
    p = p
  )
  
  # Solve the ODE system
  soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))
  
  # Subset to the last 50 years
  last_50_years <- soln[soln$time >= (max(times) - 50), ]
  
  
  # Plot the infectious population I(t) for the last 50 years
  plot(
    x = last_50_years[, "time"], 
    y = last_50_years[, "I"], 
    log = "y", 
    type = "l", 
    col = colors[i], 
    lwd = 2,
    xlab = "Time (Years)", 
    ylab = "Infectious I(t)"
  )
  abline(h = Ieqm, lty = "dotted", lwd = 2, col = "grey")
  legend(
    "topright", 
    legend = paste("p =", p), 
    col = colors[i], 
    lty = 1, 
    lwd = 2, 
    bty = "n"
  )
  
  
  # Plot the periodogram for the last 50 years
  periodogram_d(
    last_50_years, 
    xlim = c(0, 10), 
    color = colors[i]
  )
  title(main = paste("Periodogram with p =", p, "(Last 50 Years)"), line = 3)
}




##### SEE DRAMATIC CHANGES IN RANGE OF VALUES CLOSE TO 0.66 
#### ( can adapt this for what we want and choose to analyze certain ones )
## Define a finer grid of vaccination rates
p_values <- seq(0.55, 0.7, by = 0.01)
colors <- rainbow(length(p_values))

# Initialize a data frame to store results
results <- data.frame(p = numeric(), dominant_period = numeric())

# Time sequence for the simulation
times <- seq(0, 100, by = 1/52)

# Loop over vaccination rates
for (i in seq_along(p_values)) {
  p <- p_values[i]
  
  # Update parameters
  parms <- c(beta = beta, gamma = gamma, mu = mu, alpha = alpha, N = N, p = p)
  
  # Solve the ODE system
  soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))
  
  # Subset to the last 50 years
  last_50_years <- soln[soln$time >= (max(times) - 50), ]
  
  # Periodogram analysis
  s <- spectrum(last_50_years$I, plot = FALSE)
  adjusted_freq <- 52 * s$freq
  dominant_period <- 1 / adjusted_freq[which.max(s$spec)]
  
  # Store the result
  results <- rbind(results, data.frame(p = p, dominant_period = dominant_period))
  
  # Plot the time series (optional)
  plot(
    x = last_50_years$time, y = last_50_years$I, log = "y", type = "l",
    col = colors[i], lwd = 2, xlab = "Time (Years)", ylab = "Infectious I(t)",
    main = paste("Vaccination Rate p =", round(p, 2))
  )
  
  # Plot the periodogram (optional)
  plot(
    1 / (52 * s$freq), s$spec, type = "l", col = colors[i], lwd = 2,
    xlab = "Period (Years)", ylab = "Spectral Density",
    main = paste("Periodogram for p =", round(p, 2))
  )
}

# Visualize the dynamics
plot(results$p, results$dominant_period, type = "b", col = "blue", lwd = 2,
     xlab = "Vaccination Rate (p)", ylab = "Dominant Period (Years)",
     main = "Change in Dynamics with Vaccination Rate")
