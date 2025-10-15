library(LaMa)       # for HMM functions
library(RTMBdist)   # for ExGaussian distribution
library(fmesher)    # for mesh and FEM matrices
library(Matrix)     # for sparse matrices
source("./functions/forward_banded.R") # sourcing banded forward algorithm (not in LaMa yet)
source("./functions/hmm_corr.R")

### colors for plotting
color = c("black", "red", "orange")


### loading data
data <- read.csv("./data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[,c("TIME","PDCSAP_FLUX")]
colnames(data) = c("time", "y")
data <- data[2:9524,] # first observation is NA, after 9524, long series of missing

# linearly interpolating missing observations
data$y[is.na(data$y)] = approx(data$time, data$y, data$time[is.na(data$y)])$y

# centering data
data$y <- scale(data$y, scale = FALSE)

data <- data[1:4000,]

### creating mesh and finite element matrices
mesh <- fm_mesh_1d(data$time)
spde <- fm_fem(mesh)


### Simple analysis without HMM
nll <- function(par) {
  getAll(par, dat)

  ## observation model
  # parameter transformations
  sigma = exp(log_sigma); REPORT(sigma)
  f <- f0 + mu; REPORT(f0); REPORT(f) # f is the smooth function
  nll <- -sum(dnorm(y, f, sigma, log = TRUE)) # observation liklihood

  ## GP model
  # parameter transformations
  tau <- exp(log_tau); REPORT(tau)
  kappa <- exp(log_kappa); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
  omega <- plogis(logit_omega); REPORT(omega)
  Q <- tau^2 * (kappa^4 * c0 + 2 * cos(pi*omega) * kappa^2 * g1 + g2); REPORT(Q)
  nll <- nll - dgmrf(f - mu, 0, Q, log = TRUE) # GP likelihood

  nll
}

# initial parameter list
par <- list(
  log_sigma = log(7),
  log_tau = log(0.005),
  log_kappa = log(30),
  f0 = numeric(nrow(spde$c0)),
  logit_omega = qlogis(0.9),
  mu = -0.5
)

# data list
dat <- list(
  y = data$y,
  c0 = spde$c0, g1 = spde$g1, g2 = spde$g2
)

obj_simple <- MakeADFun(nll, par, random = "f0")
opt_simple <- nlminb(obj_simple$par, obj_simple$fn, obj_simple$gr)
mod_simple <- obj_simple$report()
mod_simple$tau
mod_simple$kappa
mod_simple$omega

### visualising results
idx = 1:nrow(data)
plot(data$time[idx], data$y[idx], pch = 16,
     xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Stellar flare detection")
lines(data$time[idx], mod_simple$f[idx], lwd = 2, col = "plum")
# smooth is not very periodic
# erratic behaviour whenever there is a flare


### HMM analysis
# HMM likelihood function
jnll <- function(par) {
  getAll(par, dat)

  ## state process ##
  # restricted tpm
  Gamma <- diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  # estimated initial distribution
  delta <- c(1, exp(logit_delta))
  delta <- delta / sum(delta)

  ## state-dependent process ##
  # parameter transformations
  sigma <- exp(log_sigma); REPORT(sigma)
  r <- plogis(logit_r); REPORT(r)
  lambda <- exp(log_lambda); REPORT(lambda)
  # state-dependent densities
  f <- f0 + mu; REPORT(f0); REPORT(f) # smooth function
  z <- y - f; REPORT(z); REPORT(mu)

  n <- length(z); idx = 2:n
  lallprobs <- matrix(0, n, 3)
  # regular measurement error
  lallprobs[idx,1] <- dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  lallprobs[idx,2] <- dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  lallprobs[idx,3] <- dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  # banded forward algorithm - HMM likelihood
  nll <- - forward_banded(delta, Gamma, lallprobs, k = k, logspace = TRUE)
  # lls <- forward_summands(delta, Gamma, lallprobs, logspace = TRUE)

  ### GP part ###
  # parameter transformations
  tau <- exp(log_tau); REPORT(tau)
  kappa <- exp(log_kappa); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
  omega <- plogis(logit_omega); REPORT(omega)

  Q <- tau^2 * (kappa^4 * c0 + 2 * cos(pi*omega) * kappa^2 * g1 + g2); REPORT(Q)


  nll <- nll - dgmrf(f - mu, 0, Q, log = TRUE) # GP likelihood

  nll
}

# initial parameter list
par <- list(
  eta = rep(-2, 4),
  logit_delta = rep(0, 2),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau = log(0.005),
  log_kappa = log(30),
  f0 = numeric(nrow(data)),
  logit_omega = qlogis(0.95),
  mu = -0.5
)

# data list
dat <- list(
  y = data$y,
  c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
  k = 10
)

obj <- MakeADFun(jnll, par, random = "f0")

system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)

round(p1 - opt$par, 2)

p1 <- opt$par

mod <- obj$report()
mod$kappa
mod$tau
mod$omega

states <- viterbi(mod = mod)
stateprobs <- stateprobs(mod = mod)
flare <- states != 1




sdr <- sdreport(obj, getJointPrecision = TRUE)
par <- as.list(sdr, "Est")

pacf(lls)

# cholP <- Cholesky(sdr$jointPrecision, LDL = FALSE, Imult = 1e-10)  # returns a "Cholesky" object
# pars <- lapply(1:1000, function(i){
#   z <- rnorm(nrow(sdr$jointPrecision))
#   p <- solve(cholP, z, system = "Lt")  # system="Lt" solves L^T x = z
#   p + c(sdr$par.fixed, sdr$par.random)
# })
# pars <- lapply(pars, obj3$env$parList)
#
# allstates <- matrix(NA, nrow(data), length(pars))
# allstateprobs <- array(dim = c(nrow(data), 3, length(pars)))
# for(i in 1:length(pars)){
#   getAll(pars[[i]], dat, warn = FALSE)
#   Gamma = diag(3)
#   Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
#   Gamma = Gamma / rowSums(Gamma)
#   delta = c(1, exp(logit_delta))
#   delta = delta / sum(delta)
#   sigma = exp(log_sigma)
#   r = plogis(logit_r)
#   lambda = exp(log_lambda)
#   f <- f0 + mu
#   z <- y - f
#   n <- length(z); idx = 2:n
#   allprobs = matrix(0, n, 3)
#   allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
#   allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
#   allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)
#   allstates[,i] <- viterbi(delta, Gamma, exp(allprobs))
#   allstateprobs[,,i] <- stateprobs(delta, Gamma, exp(allprobs))
# }
# state_proportions <- t(apply(allstates, 1, function(x) {
#   tab <- table(factor(x, levels = 1:3))
#   prop <- tab / length(x)
#   return(prop)
# }))
# stateprobs_unc <- apply(allstateprobs, c(1,2), mean, na.rm = TRUE)
prop_mat <- t(stateprobs)

### Plot result
# choose what to plot here
# idx <- 7000:7350
idx <- 7020:7150
# idx <- 1:nrow(data)
# prop_mat <- t(state_proportions[idx, ])  # rows = time, cols = states
# prop_mat <- t(stateprobs_unc[idx,])

# Stacked barplot
par(mfrow = c(2,1), mar = c(5,4,2,2)+0.1)
barplot(prop_mat[,idx], col = color, border = "white", space = 0,
        ylab = "State proportion",
        main = "Local state probabilities")

# Decoded time series
plot(data$time[idx], data$y[idx], col = color[states[idx]], pch = 20,
     xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Light Curve")
lines(data$time[idx], mod_simple$f[idx], lwd = 2, lty = 3, col = "blue")
lines(data$time[idx], mod$f[idx], lwd = 2, col = "plum")
legend("topleft", legend = c("Quiet", "Firing", "Decaying"), pch = 16, col = color, bty = "n")



# simulate from fitted model to check acf
set.seed(123)
Gamma <- mod$Gamma
delta <- stationary(Gamma)
sigma <- mod$sigma
r <- mod$r
lambda <- mod$lambda

nObs <- 1e5
s <- z <- rep(NA, nObs)

s[1] <- sample(1:3, 1, prob = delta)
z[1] <- 0
for(t in 2:nObs) {
  s[t] <- sample(1:3, 1, prob = Gamma[s[t-1], ])
  if(s[t] == 1){
    z[t] <- rnorm(1, 0, sigma)
  }
  if(s[t] == 2){
    z[t] <- rexgauss(1, z[t-1], sigma, lambda)
  }
  if(s[t] == 3){
    z[t] <- r * z[t-1] + rnorm(1, 0, sigma)
  }
}

acf(z)
