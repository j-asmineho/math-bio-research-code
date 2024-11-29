
### Deterministic Model
library(deSolve)

# Parameters
N <- 500000
mu <- 1 / 50                       # Birth/death rate (annual)

gamma <- 365 / 13                  # Recovery rate (13-day infectious period)
R0 <- 17                      # Basic reproduction number
beta <- R0 * (gamma + mu)    # Transmission rate
alpha <- 0.1                        # Seasonal forcing amplitude
p <- 0
constantinf<-0

# Combine into a parameter vector
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p,
  constantinf=constantinf
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
    # Seasonal forcing (alpha = 0 for no seasonality)
    beta_t <- beta * (1 + alpha * cos(2 * pi * t))
    
    # Differential equations
    dS <- mu * N - beta_t * S * I / N - mu * S - p * S
    dI <- beta_t * S * I / N - gamma * I - mu * I + constantinf
    dR <- gamma * I - mu * R + p * S - constantinf
    
    list(c(dS, dI, dR))  # Return derivatives
  })
}


# Initial conditions and time sequence
times <- seq(0, 50, by = 0.01)  
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
    s$per <- 1/(52*s$freq) # period in years (assume data are weekly)
    if (!add) { # start a new plot without plotting data
      plot(s$per,s$spec,typ="n",bty="L",ann=FALSE,las=1,
           xaxs="i",yaxs="i", xlim=xlim, col=color,...)
    }
    lines( s$per, s$spec , col=color,lwd=2)
  })
}


par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))

# Plot 
plot(
  x = soln[, "time"], 
  y = soln[, "I"], 
  log = "y", # easier to see what's going on
  type = "l", 
  col = "black", 
  lwd = 2, 
  xlab = "Time (years)", 
  ylab = "Infectious I(t)", 
  main = "Measles SIR model"
)

## add equilibrium level to see what's going on:
abline( h = Ieqm, lty = "dotted",lwd='2' , col = "grey")

legend(
  "topright", 
  legend = c("Infected (I)"), 
  col = c("black"), 
  lty = 1, 
  lwd = 2, 
  bty = "n"
)

periodogram_d(soln,trange=c(0,50),color='black')


##################################################################

p <- 0.1
constantinf<-500
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p,
  constantinf=constantinf
)
soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))

# Plot 
plot(
  x = soln[, "time"], 
  y = soln[, "I"], 
  log = "y", # easier to see what's going on
  type = "l", 
  col = "hotpink", 
  lwd = 2, 
  xlab = "Time (years)", 
  ylab = "Infectious I(t)", 
  main = paste("Measles SIR model with p = ",p)
)

## add equilibrium level to see what's going on:
abline( h = Ieqm, lty = "dotted",lwd='2' , col = "grey")

legend(
  "topright", 
  legend = c("Infected (I)"), 
  col = c("hotpink"), 
  lty = 1, 
  lwd = 2, 
  bty = "n"
)

periodogram_d(soln,trange=c(0,50),color='hotpink')
title(main=paste('Periodogram with p = ',p))


##################################################################
p <- 0.2
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p,
  constantinf=constantinf
)
soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))

# Plot 
plot(
  x = soln[, "time"], 
  y = soln[, "I"], 
  log = "y", # easier to see what's going on
  type = "l", 
  col = "dodgerblue", 
  lwd = 2, 
  xlab = "Time (years)", 
  ylab = "Infectious I(t)", 
  main = paste("Measles SIR model with p = ",p)
)

## add equilibrium level to see what's going on:
abline( h = Ieqm, lty = "dotted",lwd='2' , col = "grey")

legend(
  "topright", 
  legend = c("Infected (I)"), 
  col = c("dodgerblue"), 
  lty = 1, 
  lwd = 2, 
  bty = "n"
)

periodogram_d(soln,trange=c(0,50),color='dodgerblue')
title(main=paste('Periodogram with p = ',p))


##################################################################

p <- 0.6
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p,
  constantinf=constantinf
)
soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))

# Plot 
plot(
  x = soln[, "time"], 
  y = soln[, "I"], 
  log = "y", # easier to see what's going on
  type = "l", 
  col = "seagreen2", 
  lwd = 2, 
  xlab = "Time (years)", 
  ylab = "Infectious I(t)", 
  main = paste("Measles SIR model with p = ",p)
)

## add equilibrium level to see what's going on:
abline( h = Ieqm, lty = "dotted",lwd='2' , col = "grey")

legend(
  "topright", 
  legend = c("Infected (I)"), 
  col = c("seagreen2"), 
  lty = 1, 
  lwd = 2, 
  bty = "n"
)

periodogram_d(soln,trange=c(0,50),color='seagreen2')
title(main=paste('Periodogram with p = ',p))



#################################
p <- 0.4
parms <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  alpha = alpha,
  N = N,
  p = p,
  constantinf=constantinf
)
soln <- as.data.frame(ode(y = ic, times = times, func = SIR.vector.field, parms = parms))

# Plot 
plot(
  x = soln[, "time"], 
  y = soln[, "I"], 
  log = "y", # easier to see what's going on
  type = "l", 
  col = "purple3", 
  lwd = 2, 
  xlab = "Time (years)", 
  ylab = "Infectious I(t)", 
  main = paste("Measles SIR model with p = ",p)
)

## add equilibrium level to see what's going on:
abline( h = Ieqm, lty = "dotted",lwd='2' , col = "grey")

legend(
  "topright", 
  legend = c("Infected (I)"), 
  col = c("purple3"), 
  lty = 1, 
  lwd = 2, 
  bty = "n"
)


periodogram_d(soln,trange=c(0,50),color='purple3')
title(main=paste('Periodogram with p = ',p))

