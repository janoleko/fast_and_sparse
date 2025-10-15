hmm_cor <- function(
    Gamma, f, lag.max = 40,
    lower = -10, upper = 10, subdivisions = 1e4,
    threshold = 0.01
) {

  # stationary distribution of the Markov chain
  delta <- LaMa::stationary(Gamma)

  N <- length(delta)
  # per-state moments: m1 = E[X | state], m2 = E[X^2 | state]
  m1 <- sapply(f, function(fi) integrate(function(x) x * fi(x),
                                         lower, upper,
                                         subdivisions = subdivisions)$value)
  m2 <- sapply(f, function(fi) integrate(function(x) x^2 * fi(x),
                                         lower, upper,
                                         subdivisions = subdivisions)$value)

  # unconditional mean and second moment
  mean_X <- sum(delta * m1)
  EX2    <- sum(delta * m2)

  # variance at lag 0 (correct)
  var_X <- EX2 - mean_X^2

  # precompute Gamma^k
  Gamma_k <- vector("list", lag.max)
  Gamma_k[[1]] <- Gamma
  for (k in 2:lag.max) Gamma_k[[k]] <- Gamma_k[[k-1]] %*% Gamma

  # inner integrand for k >= 1
  inner_func <- function(x1, x2, k) {
    f1 <- sapply(f, function(fi) fi(x1))   # vector length N
    f2 <- sapply(f, function(fi) fi(x2))
    as.numeric( x1 * x2 * (((delta * f1) %*% Gamma_k[[k]]) %*% f2) )
  }

  # I1(x2,k) = integral_x1 inner_func(x1,x2,k)
  I1 <- function(x2, k) {
    integrate(function(x1) inner_func(x1, x2, k),
              lower, upper, subdivisions = subdivisions)$value
  }
  I1 <- Vectorize(I1, "x2")

  # I2(k) = double integral for lag k>=1
  # I2 <- function(k) {
  #   integrate(function(x2) I1(x2, k), lower, upper,
  #             subdivisions = subdivisions, rel.tol = 1e-10)$value
  # }
  x2seq <- seq(lower, upper, length.out = 200)
  h <- x2seq[2] - x2seq[1]
  I2 <- function(k) {
    sum(I1(xseq, k) * h)
  }

  covs <- numeric(lag.max)
  pb <- txtProgressBar(min = 1, max = lag.max, style = 3)  # style 3 = nice bar
  # lags >=1: double integral - mean^2
  for (k in 1:lag.max) {
    # cat("Computing lag", k, "\n")
    val <- I2(k)
    covs[k] <- val - mean_X^2
    setTxtProgressBar(pb, k)  # update progress bar
  }
  close(pb)

  # correlations for lags 1:lag.max
  cors <- covs / var_X

  # find lag at which correlation is smaller than threshold
  k <- which(cors < threshold)[1]

  if(is.na(k)){
    cat("Correlation does not drop below threshold. Increase lag.max\n")
  } else{
    cat("Correlation drops below", threshold, "at lag", k, "\n")
  }

  ## ---- Plot similar to acf ----
  plot(1:lag.max, cors, type = "h", lwd = 1, bty = "n",
       xlab = "Lag", ylab = "Autocorrelation",
       main = "HMM autocorrelation")
  abline(h = threshold, col = "gray")
  abline(v = k, col = "gray")
  points(k, threshold, col = "gray", pch = 20)

  invisible(list(cors = cors, delta = delta, mean = mean_X, var = var_X))
}

hmm_cor_bound <- function(Gamma, f, lag.max = 30,
                          lower = -Inf, upper = Inf, subdivisions = 1e4) {
  N <- length(f)

  # compute per-state mean and variance
  m1 <- sapply(f, function(fi)
    integrate(function(x) x * fi(x), lower, upper, subdivisions=subdivisions)$value)
  m2 <- sapply(f, function(fi)
    integrate(function(x) x^2 * fi(x), lower, upper, subdivisions=subdivisions)$value)
  sigma2 <- m2 - m1^2

  # stationary distribution
  delta <- LaMa::stationary(Gamma)

  # unconditional mean and variance
  mean_X <- sum(delta * m1)
  var_X <- sum(delta * (sigma2 + (m1 - mean_X)^2))  # total variance
  var_mu <- sum(delta * (m1 - mean_X)^2)            # variance of state means

  # second largest eigenvalue
  eig <- eigen(Gamma, symmetric = FALSE, only.values = TRUE)$values
  lambda2 <- sort(Mod(eig), decreasing = TRUE)[2]  # second largest in modulus

  # upper bound on correlations
  corr_bound <- var_mu / var_X * lambda2^(0:lag.max)

  # find lag at which correlation is smaller than threshold
  k <- which(corr_bound < threshold)[1]

  if(is.na(k)){
    cat("Correlation does not drop below threshold. Increase lag.max\n")
  } else{
    cat("Correlation drops below", threshold, "at lag", k, "\n")
  }

  # plot like acf
  plot(0:lag.max, corr_bound, type = "h", lwd = 2,
       xlab = "Lag", ylab = "Correlation bound",
       main = "HMM correlation upper bound")
  abline(h = threshold, col = "gray")
  abline(v = k, col = "gray")
  points(k, threshold, col = "gray", pch = 20)

  invisible(list(corr_bound = corr_bound, lambda2 = lambda2,
                 var_X = var_X, var_mu = var_mu, mean_X = mean_X, pi = pi))
}

