
### Deterministic Model
library(deSolve)

# Parameters
N <- 5000000
mu <- 1 / 50                       # Birth/death rate (annual)

gamma <- 365 / 13                  # Recovery rate (13-day infectious period)
R0 <- 17                      # Basic reproduction number
beta <- R0 * (gamma + mu)    # Transmission rate
alpha <- 0.1                        # Seasonal forcing amplitude
p <- NA

# Combine into a parameter vector
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p
)


## equilibrium:
R0 <- with(as.list(parms), beta/(gamma+mu))
Seqm <- (with(as.list(parms), N / R0 ))
Ieqm <- (with(as.list(parms), N  * mu/(gamma+mu)*(1 - 1/R0)))
## initial state near eqm: start near equilibrium so beginning doesn't have a huge spike
Si <- round((Seqm + (N-Seqm)*0.01))
Ii <- round((Ieqm * 0.9))
Ri <- round(N-Si-Ii)
ic <- c(S = Si, I = Ii, R = Ri)  # Initial condition vector

# SIR model with seasonal forcing and vital dynamics
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



# Initial conditions and time sequence
times <- seq(0, 50, by = 1/52)  
soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))


periodogram_d <- function(df, # data frame: date and (weekly) cases
                          xlim=c(0,10), # max 10 year period by default
                          trange=c(0,50), # whole time series by default
                          add=FALSE, color, # make new plot by default
                          ... ) {
  if (!is.na(trange[1])) { # consider only the specified time range
    df <- subset(df, time >= min(trange) & time <= max(trange))
  }
  with(df,{
    s <- spectrum(I, plot=FALSE)
    adjusted_freq <- 52*s$freq
    s$per <- 1/adjusted_freq # period in years (assume data are weekly)
    if (!add) { # start a new plot without plotting data
      plot(s$per,s$spec,typ="n",bty="L",ann=FALSE,las=1,
           xaxs="i",yaxs="i", xlim=xlim, col=color,...)
    }
    lines( s$per, s$spec , col=color,lwd=2)
  })
}


par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))


##################################################################
# Define a sequence of p values to iterate over
p_values <- c(0, 0.1, 0.4, 0.6, 0.7)

# Define colors for each p value
colors <- rainbow(length(p_values))

# Set up the plotting layout for side-by-side plots
par(mfrow = c(length(p_values), 2), mar = c(4, 4, 2, 1))  # Rows = number of p values, 2 columns

# Loop over each value of p with an index to access the corresponding color
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
  
  # Plot the infectious population I(t)
  plot(
    x = soln[, "time"],
    y = soln[, "I"],
    log = "y", # log scale to better visualize dynamics
    type = "l",
    col = colors[i],  # Use the corresponding color
    lwd = 2,
    xlab = "Time (years)",
    ylab = "Infectious I(t)",
    main = paste("Measles SIR model with p =", p)
  )
  
  # Add equilibrium line
  abline(h = Ieqm, lty = "dotted", lwd = 2, col = "grey")
  
  # Add legend
  legend(
    "topright",
    legend = paste("p =", p),
    col = colors[i],
    lty = 1,
    lwd = 2,
    bty = "n"
  )
  
  # Generate the periodogram next to the model plot
  periodogram_d(
    soln,
    trange = c(0, 50),
    color = colors[i],
    xlim = c(0, 10) # Adjust x-axis limits if needed
  )
  
  # Add title to the periodogram
  title(main = paste("Periodogram with p =", p))
}
