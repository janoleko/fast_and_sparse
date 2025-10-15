library(LaMa)
library(RTMBdist)

color = c("black", "red", "orange")

# run QFD
data <- read.csv("~/Downloads/celeriteQFD-main/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[,c("TIME","PDCSAP_FLUX")]
colnames(data) = c("time", "y")
data <- data[2:9524,] # first observation is NA, after 9524, long series of missing

# linearly interpolating missing observations
data$y[is.na(data$y)] = approx(data$time, data$y, data$time[is.na(data$y)])$y

# centering data
data$y <- scale(data$y, scale = FALSE)

pnll <- function(par) {
  getAll(par, dat)

  ## state process ##
  # restricted tpm
  Gamma = diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
  Gamma = Gamma / rowSums(Gamma)

  # estimated initial distribution
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta)

  ## state-dependent process ##
  # parameter transformations
  sigma = exp(logsigma); REPORT(sigma)
  r = plogis(logitr); REPORT(r)
  lambda = exp(loglambda); REPORT(lambda)

  # state-dependent densities
  # compute smooth
  f <- X %*% c(beta0, beta_spline); REPORT(f)
  z <- y - f; REPORT(z)

  n <- length(z)
  idx = 2:n
  allprobs = matrix(0, n, 3)
  # regular measurement error
  allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  # banded forward algorithm
  - forward(delta, Gamma, allprobs, logspace = TRUE) +
    penalty(beta_spline, S, kappa)
}

library(mgcv)
mod_gam = gam(y ~ s(time, bs = "ps", k = 200), data = data)
f_hat_gam = as.numeric(predict(mod_gam))

plot(data$time, data$y)
lines(data$time, f_hat_gam, lwd = 2, col = "orange")


# spline design matrix
k = 400
modmat = make_matrices(f ~ s(time, bs = "ps", k = k), data = data)
X = modmat$Z
S = modmat$S[[1]]

par = list(
  eta = rep(-3, 4),
  logitdelta = rep(0, 2),
  logsigma = log(7),
  logitr = qlogis(0.8),
  loglambda = log(0.02),
  beta0 = -0.5,
  beta_spline = rep(0, k-1)
)

dat = list(
  X = X,
  S = S,
  y = data$y,
  kappa = .1
)

system.time(
  mod <- qreml(pnll, par, dat,
               random = "beta_spline", psname = "kappa",
               silent = 0)
)

mod$sigma
mod$lambda
mod$r

states = viterbi(mod = mod)

plot(data$time, data$y, col = color[states], pch = 16, ylab = "y", xlab = "Time")
lines(data$time, mod$f, lwd = 3, col = "plum")

# detrended
plot(data$TIME, data$y - mod$f, col = color[states], pch = 19)


ind = 1250:1400
plot(data$TIME[ind], data$y[ind], col = color[states[ind]], pch = 16, ylim = c(-20, 40))
lines(data$TIME[ind], mod$f[ind], lwd = 4, col = "plum")


# spline modulating sine frequency

pnll2 <- function(par) {
  getAll(par, dat)

  ## state process ##
  # restricted tpm
  Gamma = diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
  Gamma = Gamma / rowSums(Gamma)

  # estimated initial distribution
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta)

  ## state-dependent process ##
  # parameter transformations
  sigma = exp(logsigma); REPORT(sigma)
  r = plogis(logitr); REPORT(r)
  lambda = exp(loglambda); REPORT(lambda)

  # state-dependent densities
  # compute smooth
  freq <- exp(X %*% c(beta0, beta_spline)); REPORT(freq)
  omega <- 2 * pi * freq; REPORT(omega)
  f <- mu + beta_mean[1] * sin(time * omega) + beta_mean[2] * cos(time * omega)
  z <- y - f; REPORT(z); REPORT(f)

  n <- length(z)
  idx = 2:n
  allprobs = matrix(0, n, 3)
  # regular measurement error
  allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  # banded forward algorithm
  - forward(delta, Gamma, allprobs, logspace = TRUE) +
    penalty(beta_spline, S, kappa)
}

# spline design matrix
k = 20
modmat = make_matrices(f ~ s(TIME, bs = "tp", k = k), data = data)
X = modmat$Z
S = modmat$S[[1]] # * 1e5

time <- data$TIME - data$TIME[1] # time in seconds
period = 0.165 # period in seconds
omega = 2 * pi / period
s <- 0.02
phi <- - omega * s
beta_mean <- 9 * c(cos(phi), sin(phi)) # mean of sine wave

plot(time, mod$f, type = "l")
lines(time, beta_mean[1] * sin(time * omega) +
        beta_mean[2] * cos(time * omega), col = "blue")

par = list(
  eta = rep(-3, 4),
  logitdelta = rep(0, 2),
  logsigma = log(7),
  logitr = qlogis(0.8),
  loglambda = log(0.02),
  beta0 = log(1/period),
  beta_spline = rep(0, k-1),
  beta_mean = beta_mean / 2,
  mu = -0.5
)

dat = list(
  X = X,
  S = S,
  y = data$y,
  time = time,
  kappa = 1e3
)

system.time(
  mod_sine <- qreml(pnll2, par, dat,
               random = "beta_spline", psname = "kappa",
               silent = 0, conv_crit = "relchange")
)

plot(data$TIME, 1/mod_sine$freq, type = "l")
abline(v = data$TIME[which(states_sine == 3)])

plot(data$TIME, data$y, pch = 16)
lines(data$TIME, mod_sine$f, lwd = 4, col = "plum")

states_sine = viterbi(mod = mod_sine)

plot(data$TIME, data$y, col = color[states], pch = 16, ylab = "y", xlab = "Time")
lines(data$TIME, mod_sine$f, lwd = 3, col = "plum")






# Using a GP for the trend ------------------------------------------------

source("./functions/forward_banded.R")

# nBasis <- 2000
# timeseq <- seq(min(data$time), max(data$time), length = nBasis)
# data <- data[1:5000,]

timeseq <- data$time
mesh <- fmesher::fm_mesh_1d(timeseq, degree = 1)
spde <- fmesher::fm_fem(mesh)
spde$X_f <- fmesher::fm_basis(mesh, data$time)

jnll <- function(par) {
  getAll(par, dat)

  ## state process ##
  # restricted tpm
  Gamma = diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
  Gamma = Gamma / rowSums(Gamma)

  # estimated initial distribution
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta)

  ## state-dependent process ##
  # parameter transformations
  sigma = exp(logsigma); REPORT(sigma)
  r = plogis(logitr); REPORT(r)
  lambda = exp(loglambda); REPORT(lambda)

  # state-dependent densities
  # compute periodic mean (parametric)
  # freq <- exp(logfreq); period <- 1 / freq; REPORT(period)
  # omega <- 2 * pi * freq
  # mu <- beta_mu[1] + beta_mu[2] * sin(time * omega) + beta_mu[3] * cos(time * omega)

  f <- f_coef + mu; REPORT(f_coef)
  # f <- as.numeric(spde$X_f %*% f_coef) + mu
  z <- y - f
  REPORT(f); REPORT(z); REPORT(mu)

  n <- length(z)
  idx = 2:n
  allprobs = matrix(0, n, 3)
  # regular measurement error
  allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  # banded forward algorithm
  nll <- - forward_banded(delta, Gamma, allprobs, k = 30, logspace = TRUE)

  #### GP part ####
  tau <- exp(logtau); REPORT(tau)
  kappa <- exp(logkappa); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1

  Q <- tau^2*(kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)        ## GMRF prior
  REPORT(Q)
  nll <- nll - dgmrf(f_coef - mu, 0, Q, log=TRUE)        ## GP likelihood

  # prior for logkappa
  # nll <- nll - dnorm(logkappa, log(3), 0.05, log = TRUE)

  nll
}

time <- data$time - data$time[1] # time in seconds
period = 0.165 # guess of period in seconds
omega = 2 * pi / period
s <- 8 # guess of shift in seconds
phi <- - omega * s # phase shift in radians
beta_mu <- 5 * c(cos(phi), sin(phi)) # corresponding coefficients

idx <- 2000:3000
plot(time[idx], data$y[idx])
lines(time[idx], beta_mu[1] * sin(time[idx] * omega) +
        beta_mu[2] * cos(time[idx] * omega), col = "blue")


par <- list(
  eta = rep(-2, 4),
  logitdelta = rep(0, 2),
  logsigma = log(7),
  logitr = qlogis(0.8),
  loglambda = log(0.025),
  logtau = log(0.015),
  logkappa = log(15),
  f_coef = numeric(nrow(data)),
  # beta_mu = c(-0.5, beta_mu),
  mu = -0.5
  # logfreq = log(1 / period)
)

dat <- list(
  y = data$y,
  spde = spde,
  time = time
)

obj <- MakeADFun(jnll, par, random = "f_coef")

# H <- obj$env$spHess(random = TRUE)
# SparseM::image(H)

system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)

mod2 <- obj$report()

# check for suitable k
eigenpairs <- eigen(mod2$Gamma)
val <- abs(eigenpairs$values[2])
cumval = cumprod(rep(val, 100))
plot(cumval, xlab = "k")
abline(v = 50)
abline(h = cumval[50])

# HMM params
mod2$sigma
mod2$r
mod2$lambda

# GP params
mod2$tau
mod2$kappa
mod2$rho


states2 <- viterbi(mod = mod2)

idx = 1:nrow(data)
# idx = 7000:9000
plot(data$time[idx], data$y[idx], col = color[states2[idx]], pch = 16,
     xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Stellar flare detection")
lines(data$time[idx], mod2$f[idx], lwd = 3, col = "plum")
legend("topright", legend = c("Quiet", "Firing", "Decaying"), pch = 16, col = color, bty = "n")

plot(data$time, mod2$mu, col = "blue", lwd = 2, type = "l", ylim = c(-11,11))
lines(data$time, mod2$f_coef, col = "lightblue",lwd = 2)
lines(data$time, mod2$mu + mod2$f_coef, col = "plum", lwd = 2)


# detrended
plot(data$time, data$y - mod2$f, col = color[states2], pch = 19)
# vs mod1
# ind <- 1700:2300
ind <- 1:nrow(data)
plot(data$time[ind], data$y[ind] - mod2$f[ind], col = color[states2[ind]], pch = 20)


# mod2
res2 <- data$y - mod2$f
idx2 <- which(states2 == 1)
hist(res2[idx2], breaks = 30, prob = TRUE)
curve(dnorm(x, sd = mod2$sigma), add = TRUE, lwd = 2, lty = 2)

# mod
res <- data$y - mod$f
idx <- which(states == 1)
hist(res[idx], breaks = 30, prob = TRUE)
curve(dnorm(x, sd = sd(res[idx])), add = TRUE, lwd = 2, lty = 2)


# sample fs
idx <- 1:500
Sigma <- as.matrix(solve(mod2$Q[idx, idx]))
mu <- mod2$mu[idx]
fs <- mvtnorm::rmvnorm(1000, mu, Sigma)







# GP with oscillating covariance ------------------------------------------

library(LaMa)
library(RTMBdist)
library(Matrix)
source("./functions/forward_banded.R")

color = c("black", "red", "orange")

# run QFD
data <- read.csv("~/Downloads/celeriteQFD-main/Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[,c("TIME","PDCSAP_FLUX")]
colnames(data) = c("time", "y")
data <- data[2:9524,] # first observation is NA, after 9524, long series of missing

# linearly interpolating missing observations
data$y[is.na(data$y)] = approx(data$time, data$y, data$time[is.na(data$y)])$y

# centering data
data$y <- scale(data$y, scale = FALSE)

# C <- spde$c0
# # Ci <- Diagonal(x = 1 / diag(C))
# G <- spde$g1
# # GCiG <- G %*% Ci %*% G
# G2 <- spde$g2
# kappa <- mod2$kappa
# tau <- mod2$tau * 10
# # omega <- 0.95
# phi <- -0.95
# # freq <- 6
# # omega <- 2 * pi * freq
# # Q <- tau^2 * (kappa^4 * C + 2 * cos(pi * omega) * kappa^2 * G + GCiG)
# Q <- tau^2 * (kappa^4 * C + 2 * phi * kappa^2 * G + G2) # GCiG)
# n <- nrow(Q)
# # Q <- tau^2 * (kappa^4 * C + 2 * kappa^2 * G + G2)        ## GMRF prior
#
# cholQ <- Cholesky(Q, LDL = FALSE, Imult = 1e-10)  # returns a "Cholesky" object
# z <- rnorm(n)
# f <- solve(cholQ, z, system = "Lt")  # system="Lt" solves L^T x = z
# plot(mod2$f, type = "l")
# lines(f)


## building FEM

mesh <- fmesher::fm_mesh_1d(data$time, degree = 1)
spde <- fmesher::fm_fem(mesh)

jnll2 <- function(par) {
  getAll(par, dat)

  ## state process ##
  # restricted tpm
  Gamma = diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
  Gamma = Gamma / rowSums(Gamma)

  # estimated initial distribution
  delta = c(1, exp(logit_delta))
  delta = delta / sum(delta)

  ## state-dependent process ##
  # parameter transformations
  sigma = exp(log_sigma); REPORT(sigma)
  r = plogis(logit_r); REPORT(r)
  lambda = exp(log_lambda); REPORT(lambda)

  # state-dependent densities
  f <- f0 + mu; REPORT(f0); REPORT(f)
  z <- y - f; REPORT(z); REPORT(mu)

  n <- length(z); idx = 2:n
  allprobs = matrix(0, n, 3)
  # regular measurement error
  allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
  # firing
  allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  # decaying
  allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)

  # banded forward algorithm
  nll <- - forward_banded(delta, Gamma, allprobs, k = 30, logspace = TRUE)

  #### GP part ####
  # parameter transformations
  tau <- exp(log_tau); REPORT(tau)
  kappa <- exp(log_kappa); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
  omega <- plogis(logit_omega); REPORT(omega)

  Q <- tau^2 * (kappa^4 * c0 + 2 * cos(pi*omega) * kappa^2 * g1 + g2); REPORT(Q)
  nll <- nll - dgmrf(f - mu, 0, Q, log = TRUE) ## GP likelihood

  nll
}

par <- list(
  eta = rep(-2, 4),
  logit_delta = rep(0, 2),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau = log(0.005),
  log_kappa = log(30),
  f0 = numeric(nrow(data)),
  logit_omega = qlogis(0.99),
  mu = -0.5
)

dat <- list(
  y = data$y,
  c0 = spde$c0, g1 = spde$g1, g2 = spde$g2
)


obj3 <- MakeADFun(jnll2, par, random = "f0")

system.time(
  opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
)

mod3 <- obj3$report()
mod3$kappa
mod3$tau
mod3$omega

states3 <- viterbi(mod = mod3)
stateprobs <- stateprobs(mod = mod3)
flare <- states3 != 1

sdr <- sdreport(obj3, getJointPrecision = TRUE)

cholP <- Cholesky(sdr$jointPrecision, LDL = FALSE, Imult = 1e-10)  # returns a "Cholesky" object
pars <- lapply(1:1000, function(i){
  z <- rnorm(nrow(sdr$jointPrecision))
  p <- solve(cholP, z, system = "Lt")  # system="Lt" solves L^T x = z
  p + c(sdr$par.fixed, sdr$par.random)
})
pars <- lapply(pars, obj3$env$parList)

allstates <- matrix(NA, nrow(data), length(pars))
allstateprobs <- array(dim = c(nrow(data), 3, length(pars)))
for(i in 1:length(pars)){
  getAll(pars[[i]], dat, warn = FALSE)
  Gamma = diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
  Gamma = Gamma / rowSums(Gamma)
  delta = c(1, exp(logit_delta))
  delta = delta / sum(delta)
  sigma = exp(log_sigma)
  r = plogis(logit_r)
  lambda = exp(log_lambda)
  f <- f0 + mu
  z <- y - f
  n <- length(z); idx = 2:n
  allprobs = matrix(0, n, 3)
  allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
  allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
  allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)
  allstates[,i] <- viterbi(delta, Gamma, exp(allprobs))
  allstateprobs[,,i] <- stateprobs(delta, Gamma, exp(allprobs))
}
state_proportions <- t(apply(allstates, 1, function(x) {
  tab <- table(factor(x, levels = 1:3))
  prop <- tab / length(x)
  return(prop)
}))
stateprobs_unc <- apply(allstateprobs, c(1,2), mean, na.rm = TRUE)

# Create stacked barplot
idx = 7000:7150
prop_mat <- t(state_proportions[idx, ])  # rows = time, cols = states
prop_mat <- t(stateprobs[idx,])
prop_mat <- t(stateprobs_unc[idx,])

# Now make stacked barplot
par(mfrow = c(2,1))
barplot(prop_mat, col = color, border = "white", space = 0,
        xlab = "Time", ylab = "State proportion",
        main = "Stacked State Proportions Over Time")

# idx = 1:nrow(data)
idx = 7000:7150
plot(data$time[idx], data$y[idx], col = color[states3[idx]], pch = 16,
     xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Stellar flare detection")
lines(data$time[idx], mod3$f[idx], lwd = 3, col = "plum")
legend("topleft", legend = c("Quiet", "Firing", "Decaying"), pch = 16, col = color, bty = "n")



cholQ <- Cholesky(mod3$Q, LDL = FALSE, Imult = 1e-10)  # returns a "Cholesky" object
z <- rnorm(nrow(mod3$Q))
f <- solve(cholQ, z, system = "Lt")  # system="Lt" solves L^T x = z

plot(data$time, f, type = "l", ylim = c(-30, 150))
lines(data$time, mod3$f, lwd = 2, col = "plum")



### model without state switching

nll <- function(par) {
  getAll(par, dat)

  # parameter transformations
  sigma = exp(log_sigma); REPORT(sigma)

  # state-dependent densities
  f <- f0 + mu; REPORT(f0); REPORT(f)

  nll <- -sum(dnorm(y, f, sigma, log = TRUE))

  #### GP part ####
  # parameter transformations
  tau <- exp(log_tau); REPORT(tau)
  kappa <- exp(log_kappa); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
  omega <- plogis(logit_omega); REPORT(omega)

  Q <- tau^2 * (kappa^4 * c0 + 2 * cos(pi*omega) * kappa^2 * g1 + g2); REPORT(Q)
  nll <- nll - dgmrf(f - mu, 0, Q, log = TRUE) ## GP likelihood

  nll
}

par <- list(
  log_sigma = log(7),
  log_tau = log(0.005),
  log_kappa = log(20),
  f0 = numeric(nrow(data)),
  logit_omega = qlogis(0.9),
  mu = -0.5
)

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

idx = 1:nrow(data)
plot(data$time[idx], data$y[idx], pch = 16,
     xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Stellar flare detection")
lines(data$time[idx], mod_simple$f[idx], lwd = 2, col = "red")
lines(data$time[idx], mod3$f[idx], lwd = 2, col = "plum")


# inferred correlation

osc_matern_corr <- function(h, kappa, theta) {
  # h: lag (can be a vector)
  # kappa: scale / damping parameter
  # theta: oscillation parameter (0 < theta < 1)

  phase <- pi * theta / 2
  damping <- kappa * cos(phase)
  oscill  <- kappa * sin(phase)

  rho <- exp(-damping * h) * sin(phase + oscill * h) / sin(phase)
  return(rho)
}

h <- seq(0, 5, length = 1000)
rho <- osc_matern_corr(h, kappa = mod3$kappa, theta = mod3$omega)
rho2 <- osc_matern_corr(h, kappa = mod_simple$kappa, theta = mod_simple$omega)

plot(h, rho, type = "l", lwd = 2, ylab = "Correlation", xlab = "Lag")
lines(h, rho2, type = "l", col = "red")


### comparing to CELERITE

celerite_sho_corr <- function(h, S0 = 1, w0 = 2*pi, Q = 2) {
  # h: lag vector
  # Defaults: S0 = 1 (amplitude scale)
  #           w0 = 2*pi (1 cycle per unit time)
  #           Q = 2 (moderate damping)
  wd <- w0 * sqrt(1 - 1/(4*Q^2))
  rho <- exp(-w0 * h / (2*Q)) * (cos(wd * h) + (1/(2*Q*wd)) * sin(wd * h))
  return(rho)
}

celerite_two_sho_corr <- function(h,
                                  S0s = c(1, 0.5),
                                  w0s = c(2*pi, 4*pi),
                                  Qs  = c(2, 5)) {
  # h: lag vector
  # Defaults:
  #   S0s = amplitudes of the two SHO components
  #   w0s = natural frequencies (rad/unit time)
  #   Qs  = quality factors (damping)

  # variance at zero lag
  k0 <- S0s[1]*w0s[1]*Qs[1]/pi + S0s[2]*w0s[2]*Qs[2]/pi

  # weighted sum of correlations
  k <- celerite_sho_corr(h, S0s[1], w0s[1], Qs[1]) * (S0s[1]*w0s[1]*Qs[1]/pi) +
    celerite_sho_corr(h, S0s[2], w0s[2], Qs[2]) * (S0s[2]*w0s[2]*Qs[2]/pi)

  rho <- k / k0
  return(rho)
}

rho2 <- celerite_two_sho_corr(h,
                            S0s = c(1, 0.5),
                            w0s = c(12.1*pi, 12.1*pi),
                            Qs  = c(100, 100))

plot(h, rho, type = "l", lwd = 2, ylab = "Correlation", xlab = "Lag")
lines(h, rho2, type = "l", col = "red")


mod3$kappa
mod3$omega
