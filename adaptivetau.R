# adaptive tau

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


library(adaptivetau)
# Transitions: each row corresponds to an event
transitions <- list(
  c(S = -1, I = +1), # Infection: S -> I
  c(I = -1, R = +1), # Recovery: I -> R
  c(S = +1),         # Birth: new susceptible
  c(S = -1),         # Death of susceptible
  c(I = -1),         # Death of infected
  c(R = -1)          # Death of recovered
)

rate_function <- function(state, parms, t) {
  with(as.list(c(state, parms)), {
    # Seasonal forcing for transmission rate
    beta_t <- beta * (1 + 0.1 * cos(2 * pi * t))
    
    rates <- c(
      beta_t * S * I / N, # Infection
      gamma * I,          # Recovery
      mu * N,             # Birth
      mu * S,             # Death of susceptible
      mu * I,             # Death of infected
      mu * R              # Death of recovered
    )
    return(rates)
  })
}
# Time limits
tmax <- 50  # Maximum time in years

# Run simulation
result <- ssa.adaptivetau(
  init.values = ic,
  transitions = transitions,
  rateFunc = rate_function,
  params = parms,
  tf = tmax
)
# Convert result to data frame for plotting
result_df <- as.data.frame(result)
colnames(result_df) <- c("time", "S", "I", "R")


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

# Plot infectious individuals over time
points(
  result_df$time, result_df$I,
  type = "l", col = "black", lwd = 1,
  xlab = "Time (years)", ylab = "Infectious I(t)",
  main = "SIR Model Simulation with Adaptive Tau"
)
