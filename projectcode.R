### Deterministic Model
library(deSolve)

# Parameters
N <- 5000000
mu <- 1 / 50                       # Birth/death rate (annual)

gamma <- 365 / 13                  # Recovery rate (13-day infectious period)
R0 <- 17                      # Basic reproduction number
beta <- R0 * (gamma + mu)    # Transmission rate
alpha <- 0.1                         # Seasonal forcing amplitude

# Combine into a parameter vector
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N
)

## equilibrium:
R0 <- with(as.list(parms), beta/(gamma+mu))
R0
Seqm <- (with(as.list(parms), N / R0 ))

Ieqm <- (with(as.list(parms), N  * mu/(gamma+mu)*(1 - 1/R0)))
## initial state near eqm: start near equilibrium so beginning doesn't have a huge spike
Si <- round((Seqm + (N-Seqm)*0.01))
Ii <- round((Ieqm * 0.9))
Ri <- round(N-Si-Ii)
ic <- c(S = Si, I = Ii, R = Ri)  # Initial condition vector

ic

# SIR model with seasonal forcing and vital dynamics
SIR.vector.field <- function(t, vars, parms) {
  with(as.list(c(parms, vars)), {
    # Seasonal forcing (alpha = 0 for no seasonality)
    beta_t <- beta * (1 + alpha * cos(2 * pi * t))
    
    # Differential equations
    dS <- mu * N - beta_t * S * I / N - mu * S
    dI <- beta_t * S * I / N - gamma * I - mu * I
    dR <- gamma * I - mu * R
    
    list(c(dS, dI, dR))  # Return derivatives
  })
}


# Initial conditions and time sequence
times <- seq(0, 50, by = 0.01)  
soln <- ode(y = ic, times = times, func = SIR.vector.field, parms = parms)

# Plot 
plot(
  x = soln[, "time"], 
  y = soln[, "I"], 
  log = "y", # easier to see what's going on
  type = "l", 
  col = "red", 
  lwd = 2, 
  xlab = "Time (years)", 
  ylab = "Infectious I(t)", 
  main = "Measles SIR model"
)

## add equilibrium level to see what's going on:
abline( h = Ieqm, lty = "dotted", col = "grey")

legend(
  "topright", 
  legend = c("Infected (I)"), 
  col = c("red"), 
  lty = 1, 
  lwd = 2, 
  bty = "n"
)


### Stochastic Model
SIR.Gillespie <- function(parms, ic, tmax= 50, dtsave = 1/52, yearstep = 1) {
  start.time <- proc.time()
  t <- 0 # start at time 0
  tvec <- c(t) # vector of event times
  beta = parms['beta']
  mu = parms['mu']
  gamma = parms['gamma']
  N=parms['N']
  
  # equilibrium points (initialized above)
  # Sstar <- Si
  # Istar <- Ii
  
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
    newinfected <- 0   # Add new infected to add noise
    
    # Total event rate possible
    a0 <- infectionrate + recoveryrate + birthrate + deathrateS + deathrateI + deathrateR + newinfected
    
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
      # R <- R + 1 # just to keep same population
    } else if (r < infectionrate + recoveryrate + birthrate + deathrateS + deathrateI) { # Death of an infected occurs
      I <- I - 1
      # R <- R + 1
    } else if (r < infectionrate + recoveryrate + birthrate + deathrateS + deathrateI + deathrateR) {
      R <- R - 1 # no change
    } else {
      I <- I + newinfected
      R <- R - newinfected
    }
    
    if (t >= tsave){
      # Store updated results
      tvec <- c(tvec, t)
      Svec <- c(Svec, S)
      Ivec <- c(Ivec, I)
      Rvec <- c(Rvec, R)
      
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
  print(data.frame(S=S,I=I,R=R))
  end.time <- proc.time() - start.time
  message("run time: ")
  print(end.time)
  
  # Create data frame of results
  df<-data.frame(time = tvec, S = Svec, I = Ivec, R = Rvec)
  if (a0==0) df <- df[-nrow(df),]
  cat(nrow(df), " events returned by SI.Gillespie.\n")
  return(df)
}

# Function to plot stochastic simulations and the periodogram
plot_multiple_gillespie_lines <- function(parms, ic, tmax = 50, nsim = 5) {
  start.time <- proc.time()
  
  # Run multiple simulations
  result_list <- lapply(1:nsim, function(x) SIR.Gillespie(parms = parms, ic = ic, tmax = tmax))
  
  # Determine plot limits based on all simulations
  y_range <- range(sapply(result_list, function(res) res$I), na.rm = TRUE)
  
  # Create base plot
  plot(NULL, xlim = c(0, tmax), ylim = y_range,
       xlab = "Time (years)", ylab = "Infectious I(t)", 
       main = "Stochastic Gillespie Simulations")
  
  colors <- rainbow(nsim)
  
  # Plot each simulation with its unique color
  for (i in 1:nsim) {
    lines(result_list[[i]]$time, result_list[[i]]$I, col = colors[i], lwd = 1.5)
    print(paste("Simulation", i, "complete."))
  }
  
  # plot deterministic model (if available)
  lines(x = soln[, "time"], y = soln[, "I"], col = "black", lwd = 2)
  
  # Add legend
  legend("topright", legend = paste("Simulation", 1:nsim), 
         col = colors, lty = 1, bty = "n", cex = 0.7, lwd = 1.5)
  
  # Now, calculate and plot the periodogram for one of the simulations (e.g., the first one)
  v <- result_list[[1]]$I  # Use the 'I' time series from the first simulation
  
  # Calculate the periodogram using the 'spectrum' function
  s <- spectrum(v, plot = FALSE)
  
  # Adjust the frequency to be in terms of years (optional)
  adjusted_freq <- s$freq * 52  # if data is weekly, convert to yearly frequencies
  
  # Create the periodogram plot
  plot(1 / adjusted_freq, s$spec, type = "l", col = "blue", xlim= c(0,5),
       xlab = "Years", ylab = "Periodogram", main = "Periodogram of Infectious Population")
  
  end.time <- proc.time() - start.time
  message("run time: ")
  print(end.time)
}

# Run the function to plot the stochastic simulations and periodogram
plot_multiple_gillespie_lines(parms = parms, ic = ic)





