
## Total run time: ~ 5.5 hours (can reduce run time in stochastic model below by selecting less vaccination rates (p))


##### Allow to settle after times 1000, then plot only for the last 50 years
library(deSolve)

# Parameters
N <- 5000000                      # Total population
mu <- 1 / 50                      # Birth/death rate (annual)
gamma <- 365 / 13                 # Recovery rate (13-day infectious period)
R0 <- 17                          # Basic reproduction number
beta <- R0 * (gamma + mu)         # Transmission rate
alpha <- 0.1                      # Seasonal forcing amplitude
p <- 0                           # Placeholder for vaccination rate

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
periodogram <- function(df, xlim = c(0, 10), color, ...) {
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
  periodogram(
    last_50_years, 
    xlim = c(0, 10), 
    color = colors[i]
  )
  title(main = paste("Periodogram with p =", p, "(Last 50 Years)"), line = 3)
}




##### SEE DRAMATIC CHANGES IN RANGE OF VALUES CLOSE TO 0.66 
#### ( can adapt this for what we want and choose to analyze certain ones )
## Define a finer grid of vaccination rates
# Define a finer grid of vaccination rates
p_values <- seq(0, 1, by = 0.01)
# colors <- rainbow(length(p_values))

# Initialize a data frame to store results
results <- data.frame(p = numeric(), dominant_period = numeric())

# Time sequence for the simulation
times <- seq(0, 200, by = 1/52)

# Open a PDF device to save the plots
# pdf("/Users/sarah/Desktop/McMaster/Math4MB3/project/images/Changes in Vaccination Rates p = 0 to 1.pdf", width = 8, height = 4)

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
  
  # Set up side-by-side plots
  par(mfrow = c(1, 2))
  
  # Plot the time series dynamics (left graph)
  plot(
    x = last_50_years$time, y = last_50_years$I, type = "l",
    col = "red", lwd = 2, xlab = "Time (Years)", ylab = "Infectious I(t)",
    main = paste("Dynamics for p =", round(p, 2))
  )
  
  # Plot the periodogram (right graph)
  plot(
    1 / (52 * s$freq), s$spec, type = "l", col = "red", lwd = 2,
    xlab = "Period (Years)", ylab = "Spectral Density",
    xlim = c(0, 6),
    main = paste("Periodogram for p =", round(p, 2))
  )
}

par(mfrow = c(1, 1))
# Visualize the dynamics
plot(
  results$p, results$dominant_period, type = "b", col = "blue", lwd = 2, pch = 20,
  xlab = "Vaccination Rate (p)", ylab = "Dominant Period (Years)",
  main = "Change in Dynamics with Vaccination Rate (Deterministic Model)"
)

# Close the PDF device
# dev.off()





# Stochastic model, start at end values from deterministic model once settled 
### 15.5 minutes
# reset parms
p <- 0.66

# Combine into a parameter vector
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p
)

soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))

# set initial value of stochastic model to final value from deterministic model
final_state <- tail(soln, 1)
ic_stochastic <- c(S = round(final_state$S), I = round(final_state$I), R = round(final_state$R))

# Number of simulations for Gillespie simulation
nsim <- 5
tmax <- 100


SIR.Gillespie <- function(parms, ic, tmax, dtsave = 1/52, yearstep = 1) {
  start.time <- proc.time()
  t <- 0 # start at time 0
  tvec <- c(t) # vector of event times
  beta = parms['beta']
  mu = parms['mu']
  gamma = parms['gamma']
  N=parms['N']
  p=as.numeric(parms['p'])
  
  # Begin closer to eqm points
  S <- ic['S'] 
  I <- ic['I']
  R <- ic['R']
  
  message("Initial state: S I R")
  print(data.frame(S=S,I=I,R=R))
  
  # initialize SIR vectors
  Svec <- c(S)              
  Ivec <- c(I)               
  Rvec <- c(R)
  
  # set the time we save as the first time step dtsave
  tsave <- dtsave
  
  # set year output when printing
  year <- yearstep
  
  
  while (t <= tmax) {
    
    # Compute seasonal transmission rate
    beta_t <- beta * (1 + alpha * cos(2 * pi * t))
    
    # Compute event rates
    infectionrate <- beta_t * S * I / N       # Rate of infection
    recoveryrate <- gamma * I              # Rate of recovery
    birthrate <- mu * N                    # Rate of births
    deathrateS <- mu * S                  # Death rate for susceptibles
    deathrateI <- mu * I                  # Death rate for infected
    deathrateR <- mu * R                  # Death rate for recovered
    newbornvaccinatedrate <- mu * p * N       # Rate of vaccinated newborns
    newinfected <- 0   # Add new infected to add noise
    
    # Total event rate possible
    a0 <- infectionrate + recoveryrate + birthrate + deathrateS + deathrateI + deathrateR + newinfected + newbornvaccinatedrate
    
    ## compute time to next event
    dt <- (1/a0)*log(1/runif(1))
    t <- t + dt
    
    # find random value
    r <- runif(1,min=0,max=a0)
    # see which rate occurs
    if (r < infectionrate) { # Infection occurs
      S <- S - 1
      I <- I + 1
    } else if (r < infectionrate + recoveryrate) { # Recovery occurs
      I <- I - 1
      R <- R + 1
    } else if (r < infectionrate + recoveryrate + birthrate) { # Birth occurs
      S <- S + 1
      R <- R - 1 # just to keep same population
    } else if (r < infectionrate + recoveryrate + birthrate + deathrateS) { # Death of a susceptible occurs
      S <- S - 1
      R <- R + 1 # just to keep same population
    } else if (r < infectionrate + recoveryrate + birthrate + deathrateS + deathrateI) { # Death of an infected occurs
      I <- I - 1
      R <- R + 1
    } else if (r < infectionrate + recoveryrate + birthrate + deathrateS + deathrateI + deathrateR) {
      R <- R # no change
    } else if (r < infectionrate + recoveryrate + birthrate + deathrateS + deathrateI + deathrateR + newinfected){
      I <- I + newinfected
      R <- R - newinfected
    } else {
      S <- S - 1
      R <- R + 1
    }
    
    if (t >= tsave){
      # Store updated results
      tvec <- c(tvec, t)
      Svec <- c(Svec, S)
      Ivec <- c(Ivec, I)
      Rvec <- c(Rvec, R)
      # print(data.frame(S=S,I=I,R=R))
      
      tsave <- tsave + dtsave
      
      if (t >= year){
        # Print each year when running to see how quickly it runs
        message('year')
        print(tsave)
        year <- yearstep + year
      }
    }
  }
  
  message("final state: ")
  end.time <- proc.time() - start.time
  message("run time: ")
  print(end.time)
  
  # Create data frame of results
  df<-data.frame(time = tvec, S = Svec, I = Ivec, R = Rvec)
  if (a0==0) df <- df[-nrow(df),]
  cat(nrow(df), " events returned by SI.Gillespie.\n")
  return(df)
}

### Function to plot stochastic simulations and the periodogram
# Takes approximately 10 minutes to run with N = 500,000 and tmax = 50

par(mfrow = c(1, 1))
plot_multiple_gillespie_lines <- function(parms, ic, tmax, nsim) {
  start.time <- proc.time()
  
  p <- as.numeric(parms['p'])
  
  # Run multiple simulations and print simulation complete to track progress
  result_list <- lapply(1:nsim, function(x) {
    set.seed(x)  # Set a unique seed for each simulation based on its index
    res <- SIR.Gillespie(parms = parms, ic = ic_stochastic, tmax = tmax)
    message(sprintf("Simulation %d completed for p = %.2f", x, p))
    return(res)
  })
  
  
  last_50_start_time <- tmax - 50
    
  last_50_deterministic <- soln[soln$time >= last_50_start_time, ]
    
    # Determine plot y limits based on all simulations
  y_range <- range(c(
      sapply(result_list, function(res) res$I),
      last_50_deterministic$I
    ), na.rm = TRUE)
    
    # Create base plot
  plot(NULL, xlim = c(tmax-50, tmax), ylim = y_range,
       xlab = "Time (years)", ylab = "Infectious I(t)", 
       main = paste("Stochastic and Deterministic Models (last 50 years)\n", "Vaccination Rate (p) = ", p))
  
  colors <- rainbow(nsim)
  
  simulation_data <- data.frame(simulation = integer(0), time = numeric(0), I = numeric(0))
  
  for (i in 1:nsim) {
    # Subset the stochastic simulation to the last 50 years
    last_50_stochastic <- result_list[[i]][result_list[[i]]$time >= last_50_start_time, ]
    
    # Plot the filtered stochastic results
    lines(last_50_stochastic$time, last_50_stochastic$I, col = colors[i], lwd = 1.5)
    
    simulation_data <- rbind(simulation_data, 
                             data.frame(simulation = rep(i, nrow(last_50_stochastic)), 
                                        time = last_50_stochastic$time, 
                                        I = last_50_stochastic$I))
  }
    
    
  # plot deterministic model (if available)
  lines(x = last_50_deterministic[, "time"], y = last_50_deterministic[, "I"], col = "black", lwd = 2)
  
  # Add legend
  legend("topright", legend = paste("Simulation", 1:nsim), 
         col = colors, lty = 1, bty = "n", cex = 0.7, lwd = 1.5)
  
  # Calculate and plot the periodogram for each simulation
  par(mfrow = c(3, 2))  # Arrange plots in a grid (5 stochastic + 1 deterministic)
  
  # Loop through the simulations
  for (i in 1:nsim) {
    # Subset the data for the current simulation (simulation 'i')
    sim_data <- subset(simulation_data, simulation == i)
    
    # Extract the 'I' values from the simulation data
    v <- sim_data$I
    
    # Calculate the periodogram using the 'spectrum' function
    s <- spectrum(v, plot = FALSE)
    
    # Adjust the frequency to be in terms of years (optional, based on the frequency of your data)
    adjusted_freq <- s$freq * 52  # If data is weekly, convert to yearly frequencies
    
    # Create the periodogram plot for the current simulation
    plot(1 / adjusted_freq, s$spec, type = "l", col = colors[i],
         xlim = c(0, 5), xlab = "Years", ylab = "Periodogram",
         main = sprintf("Periodogram of Simulation %d", i))
  }
  
  # Deterministic periodogram (bottom right)
  deterministic_s <- spectrum(last_50_deterministic$I, plot = FALSE)
  adjusted_freq_d <- deterministic_s$freq * 52  # Convert to yearly frequencies
  plot(1 / adjusted_freq_d, deterministic_s$spec, type = "l", col = "black",
       xlim = c(0, 5), xlab = "Years", ylab = "Periodogram",
       main = sprintf("Deterministic Periodogram for p = %.2f", p))


  par(mfrow = c(1, 1))  # Reset plot layout to single panel
  
  end.time <- proc.time() - start.time
  message("total run time: ")
  print(end.time)
  
}



# Run the function to plot the stochastic simulations and periodogram
### Can run to see for p = 0.66 (annual cycles)
plot_multiple_gillespie_lines(parms = parms, ic = ic, tmax = tmax, nsim = nsim)

# Define a sequence of p values
p_values <- c(0.1, 0.11, 0.35, 0.36, 0.45, 0.46, 0.57, 0.58, 0.6, 0.61, 0.67, 0.68, 0.69, 0.72, 0.73, 0.74, 0.87, 0.91, 0.93)  # Adjust the step size for desired resolution

# Time sequence for the simulation
times <- seq(0, 100, by = 1/52)


#### NOTE: takes about 5.5 hours to load, can decrease number of p_values for less run time

start.time <- proc.time()

# Loop over p values
for (p in p_values) {
  # Update the parameters for the current p value
  parms <- c(
    beta = beta,
    gamma = gamma,
    mu = mu,
    alpha = alpha,
    N = N,
    p = p
  )
  
  soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))
  
  # set initial value of stochastic model to final value from deterministic model
  final_state <- tail(soln, 1)
  ic_stochastic <- c(S = round(final_state$S), I = round(final_state$I), R = round(final_state$R))
  
  # Save the output to a PDF for each p value
  pdf(file = paste0("/Users/sarah/Desktop/McMaster/Math4MB3/project/images/stochasticsim_vaccinationrate_", formatC(p, digits = 2, format = "f"), ".pdf"))
  
  # Run the function
  plot_multiple_gillespie_lines(parms = parms, ic = ic, tmax = tmax, nsim = nsim)
  
  dev.off()  # Close the PDF device
}
end.time <- proc.time() - start.time
message("total run time: ")
print(end.time)
