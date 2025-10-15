library(LaMa) # also loads RTMB

## compute marginal state distribution at each time point
statedist <- function(delta, Gamma) {
  nObs <- dim(Gamma)[3] # number of observations
  nStates <- dim(Gamma)[1] # number of states

  # Initialise Delta matrix to store state distributions
  Delta <- matrix(NaN, nrow = nObs, ncol = nStates)

  # Set initial state distribution
  Delta[1, ] <- as.numeric(delta)

  # Loop through each observation
  for (t in 2:nObs) {
    # Compute the state distribution at time t
    Delta[t, ] <- Delta[t - 1, ] %*% Gamma[, , t]
    # Delta[t, ] <- t(Gamma[, , t]) %sp% Delta[t - 1, ]
  }

  return(Delta)
}

# about twice the cost - but banded if Gamma depends on random effect
statedist_banded <- function(delta, Gamma, j = 30){
  nObs <- dim(Gamma)[3] # number of observations
  nStates <- dim(Gamma)[1] # number of states

  n_tracks = ceiling(nObs / j) - 2 # number of remaining tracks
  last_start = n_tracks * j + j + 1

  # fixed intermediate state distribution
  Rho <- rep(1, nStates) / nStates

  # Initialise Delta matrix to store state distributions
  Delta <- matrix(NaN, nrow = nObs, ncol = nStates)

  # Set initial state distribution
  Delta[1, ] <- as.numeric(delta)

  for(t in 2:(2*j)){
    # Compute the state distribution at time t
    Delta[t, ] <- Delta[t - 1, ] %*% Gamma[, , t]
    # Delta[t, ] <- t(Gamma[, , t]) %sp% Delta[t - 1, ]
  }

  startInd <- j + 1
  endInd <- j + 2 * j

  for(track in seq_len(n_tracks)){
    # foo <- Rho
    foo <- Rho

    end <- min(endInd, nObs)
    for(t in startInd:end){
      foo <- foo %*% Gamma[, , t]
      # foo <- foo %sp% Gamma[, , t]
      if(t > startInd + j - 1){
        Delta[t, ] <- foo
      }
    }
    startInd <- startInd + j
    endInd <- endInd + j
  }
  Delta
}

# banded version of the transition probability matrix
# this is basically the exact same as in LaMa,
# with the only change that it uses %sp% instead of %*% because the latter
# currently kills automatic sparsity detection for some reason
tpm_g_banded = function(Z, beta, byrow = FALSE, report = TRUE){

  K = nrow(beta)
  p = ncol(beta) - 1
  N = as.integer(0.5 + sqrt(0.25 + K), 0)

  Z = as.matrix(Z)

  if(ncol(Z) == p){
    Z = cbind(1, Z) # adding intercept column
  } else if(ncol(Z) != p + 1){
    stop("The dimensions of Z and beta do not match.")
  }

  # report quantities for easy use later
  if(report) {
    # Setting colnames for beta: Inherit colnames from Z
    colnames(beta) <- colnames(Z)
    if(is.null(rownames(beta))){
      # Setting rownames: depends on byrow
      names <- outer(paste0("S", 1:N, ">"), paste0("S", 1:N), FUN = paste0) # matrix
      diag(names) <- NA
      rownames(beta) <- na.omit(if (byrow) c(t(names)) else c(names))
    }
    RTMB::REPORT(beta)
  }

  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")

  expEta = exp(Z %sp% t(beta))
  Gamma = array(1, dim = c(N, N, nrow(expEta)))

  ## Loop over entries (stuff over time happens vectorised which speeds up the tape)
  col_ind <- 1
  for(i in seq_len(N)){ # loop over rows
    for(j in seq_len(N)){ # loop over columns
      if(j != i){ # only if non-diagonal
        if(byrow){
          Gamma[i, j, ] <- expEta[, col_ind]
        } else{
          Gamma[j, i, ] <- expEta[, col_ind]
        }
        # increase col_ind by one
        col_ind = col_ind + 1
      }
    }
  }

  # Normalise rows to sum to 1
  for(i in seq_len(N)){
    # transposing is necessary because Gamma[i,,] / colSums(Gamma[i,,]) does not work as expected
    Gamma[i, , ] <- t(t(Gamma[i, , ]) / rowSums(t(Gamma[i, , ])))
  }

  # naming
  statenames <- paste0("S", 1:N)
  rownames(Gamma) <- statenames
  colnames(Gamma) <- statenames
  Gamma
}

# banded version of the continuous-time transition probability matrix
tpmCT_g_banded <- function(Z, beta, dt, byrow = FALSE, report = TRUE){
  K <- nrow(beta)
  p <- ncol(beta) - 1
  N <- as.integer(0.5 + sqrt(0.25 + K))

  Z <- as.matrix(Z)
  if (ncol(Z) == p) Z <- cbind(1, Z)
  if (ncol(Z) != p + 1) stop("The dimensions of Z and beta do not match.")

  if (report) {
    colnames(beta) <- colnames(Z)
    if (is.null(rownames(beta))) {
      names <- outer(paste0("S", 1:N, ">"), paste0("S", 1:N), FUN = paste0)
      diag(names) <- NA
      rownames(beta) <- na.omit(if (byrow) c(t(names)) else c(names))
    }
    RTMB::REPORT(beta)
  }

  "[<-" <- ADoverload("[<-")
  "c"   <- ADoverload("c")

  # expEta: T x K (K = N*(N-1))
  expEta <- exp(Z %sp% t(beta))
  Tlen   <- nrow(expEta)

  Gamma <- array(0, dim = c(N, N, Tlen))
  diag(Gamma[, , 1]) <- 1

  # Fast for N=2 using closed form
  if (N == 2L) {

    tiny <- 1e-12
    for (t in 2:Tlen) {

      a <- expEta[t, 1]
      b <- expEta[t, 2]
      s <- a + b
      h <- dt[t]

      s_eps <- s + tiny
      e <- exp(-s * h)

      P11 <- (b + a * e) / s_eps
      P12 <- (a - a * e) / s_eps
      P21 <- (b - b * e) / s_eps
      P22 <- (a + b * e) / s_eps

      Gamma[1,1,t] <- P11; Gamma[1,2,t] <- P12
      Gamma[2,1,t] <- P21; Gamma[2,2,t] <- P22
    }
    dimnames(Gamma) <- list(paste0("S",1:N), paste0("S",1:N), NULL)
    return(Gamma)
  }

  # Generic N fallback – tape-safe scalar writes, no AD comparisons
  Q <- matrix(0, N, N)

  # Precompute k -> (i,j) pairs
  ij_i <- integer(K)
  ij_j <- integer(K)
  k <- 1L
  if (byrow) {
    for (i in 1:N) for (j in 1:N) if (j != i) { ij_i[k] <- i; ij_j[k] <- j; k <- k + 1L }
  } else {
    for (j in 1:N) for (i in 1:N) if (i != j) { ij_i[k] <- i; ij_j[k] <- j; k <- k + 1L }
  }

  for (t in 2:Tlen) {
    Q[,] <- 0

    if (byrow) {
      for (k in 1:K) {
        i <- ij_i[k]; j <- ij_j[k]
        val <- expEta[t, k]
        Q[i, j] <- val
        Q[i, i] <- Q[i, i] - val
      }
    } else {
      for (k in 1:K) {
        i <- ij_i[k]; j <- ij_j[k]
        val <- expEta[t, k]
        Q[j, i] <- val
        Q[i, i] <- Q[i, i] - val
      }
    }
    Gamma[,,t] <- as.matrix(Matrix::expm(Q * dt[t]))
  }

  dimnames(Gamma) <- list(paste0("S",1:N), paste0("S",1:N), NULL)
  Gamma
}


## the approximation here works like
# f(x_1, ..., x_t) =
# f(x_1, ..., x_k) * f(x_k+1, ..., x_2k | x_1, ..., x_k) * f(x_2k+1, ..., x_3k | x_1, ..., x_2k) *
# f(x_3k+1, ..., x_4k | x_1, ..., x_3k) * ... (true likelihood)
# we truncate the conditioning to get:
# f(x_1, ..., x_t) ~=
# f(x_1, ..., x_k) * f(x_k+1, ..., x_2k | x_1, ..., x_k) * f(x_2k+1, ..., x_3k | x_k+1, ..., x_2k) *
# f(x_3k+1, ..., x_4k | x_2k+1, ..., x_3k) * ... (approximate likelihood)
# if k is chosen sufficiently large, there is basically no difference

## banded version of the forward algorithm with one tpm
forward_banded <- function(delta, # initial distribution
                           Gamma, # tpm array
                           allprobs, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                           k = 20, # truncation parameter for sparse likelihood approximation
                           report = TRUE, # report quantities of interest?
                           logspace = FALSE # are allprobs on log-scale? -> internal computations in logspace
                           ){
  # AD overloading (currently necessary with RTMB)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  if(logspace){
    logallprobs <- allprobs # assuming allprobs is on log-scale
    allprobs <- exp(logallprobs) # to report on the original scale
  }

  if(report){
    # reporting
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
  }

  # number of observations
  nObs <- nrow(allprobs)

  # assigning initial distribution as delta1
  # making delta1 a matrix for compatibility with RTMB
  delta1 <- matrix(delta, nrow = 1, ncol = length(delta))

  # overwriting delta as stationary
  delta <- stationary(Gamma)
  # making delta a matrix for compatibility with RTMB
  delta <- matrix(delta, nrow = 1, ncol = length(delta))

  if(!logspace){
    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    foo <- delta1 * allprobs[1, , drop = FALSE]
    sumfoo <- sum(foo)
    l <- log(sumfoo)
    phi <- foo / sumfoo
    for(t in 2:(2*k)){
      foo <- (phi %*% Gamma) * allprobs[t, , drop = FALSE]
      # foo <- (phi %sp% Gamma) * allprobs[t, , drop = FALSE]
      sumfoo <- sum(foo)
      l <- l + log(sumfoo)
      phi <- foo / sumfoo
    }

    # now we always compute f(x_t, ... x_t+k | x_t-1, ..., x_t-k)
    # = f(x_t-k, ..., x_t-1, x_t, ... x_t+k) / f(x_t-k, ..., x_t-1)
    # the nice thing is that both the numerator and denominator can be
    # computed with one pass of the forward algorithm
    startInd <- k + 1
    endInd <- k + 2*k

    n_tracks = ceiling(nObs / k) - 2 # number of remaining tracks
    for(track in seq_len(n_tracks)){
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      foo <- delta * allprobs[startInd, , drop = FALSE]
      sumfoo <- sum(foo)
      log_num <- log(sumfoo)
      phi <- foo / sumfoo
      end <- min(endInd, nObs)
      for(t in (startInd + 1):end){
        foo <- (phi %*% Gamma) * allprobs[t, , drop = FALSE]
        # foo <- (phi %sp% Gamma) * allprobs[t, , drop = FALSE]
        sumfoo <- sum(foo)
        log_num <- log_num + log(sumfoo)
        phi <- foo / sumfoo
        # at the middle index, save value for denominator
        if(t == (startInd + k - 1)){
          log_denom <- log_num
        }
      }
      l <- l + log_num - log_denom

      startInd <- startInd + k
      endInd <- endInd + k
    }
  } else {
    delta1 <- delta1 + 1e-8
    delta1 <- delta1 / sum(delta1) # to avoid numerical issues

    delta <- delta + 1e-8
    delta <- delta / sum(delta) # to avoid numerical issues

    logdelta1 <- matrix(log(delta1), nrow = 1)
    logdelta <- matrix(log(delta), nrow = 1)

    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    logfoo <- logdelta1 + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:(2*k)){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
      # logfoo <- log(exp(logfoo) %sp% Gamma) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l <- l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }

    # now we always compute f(x_t, ... x_t+k | x_t-1, ..., x_t-k)
    # = f(x_t-k, ..., x_t-1, x_t, ... x_t+k) / f(x_t-k, ..., x_t-1)
    # the nice thing is that both the numerator and denominator can be
    # computed with one pass of the forward algorithm
    startInd <- k + 1
    endInd <- k + 2*k

    n_tracks = ceiling(nObs / k) - 2 # number of remaining tracks
    for(track in seq_len(n_tracks)){
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      logfoo <- logdelta + logallprobs[startInd, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      log_num <- logsumfoo
      logfoo <- logfoo - logsumfoo
      end <- min(endInd, nObs)
      for(t in (startInd + 1):end){
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        # logfoo <- log(exp(logfoo) %sp% Gamma) + logallprobs[t, , drop = FALSE]
        logsumfoo <- LaMa:::logspace_add(logfoo)
        log_num <- log_num + logsumfoo
        logfoo <- logfoo - logsumfoo
        # at the middle index, save value for denominator
        if(t == (startInd + k - 1)){
          log_denom <- log_num
        }
      }
      l <- l + log_num - log_denom

      startInd <- startInd + k
      endInd <- endInd + k
    }
  }

  # returning approximate log-likelihood
  l
}


forward_banded2 <- function(delta, # initial distribution
                            Gamma, # tpm array
                            allprobs, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                            k = 20, # truncation parameter for sparse likelihood approximation
                            report = TRUE, # report quantities of interest?
                            logspace = FALSE # are allprobs on log-scale? -> internal computations in logspace
){
  # AD overloading (currently necessary with RTMB)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  # TapeConfig(matmul = "compact")

  if(logspace){
    logallprobs <- allprobs # assuming allprobs is on log-scale
    allprobs <- exp(logallprobs) # to report on the original scale
  }

  if(report){
    # reporting
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
  }

  if(length(k) == 1) k <- rep(k, 2)
  if(k[2] > k[1]) stop("k[2] must be smaller than or equal to k[1].")
  if(k[2] < 1) stop("k[2] must be at least 1.")

  # number of observations
  nObs <- nrow(allprobs)

  # assigning initial distribution as delta1
  # making delta1 a matrix for compatibility with RTMB
  delta1 <- matrix(delta, nrow = 1, ncol = length(delta))
  # overwriting delta as stationary
  delta <- stationary(Gamma)
  # making delta a matrix for compatibility with RTMB
  delta <- matrix(delta, nrow = 1, ncol = length(delta))

  if(!logspace){
    # initialise: compute f(x_1, ..., x_(k[1] + k[2]))
    # regular forward algorithm:
    foo <- delta1 * allprobs[1, , drop = FALSE]
    sumfoo <- sum(foo)
    l <- log(sumfoo)
    foo <- foo / sumfoo
    for(t in 2:sum(k)){
      foo <- (foo %*% Gamma) * allprobs[t, , drop = FALSE]
      sumfoo <- sum(foo)
      l <- l + log(sumfoo)
      foo <- foo / sumfoo
    }

    # now loop over tracks
    current_t <- sum(k) # t
    startInd <- current_t - k[1] + 1 # t - k[1] + 1
    endInd <- current_t + k[2] # t + k[2]

    while(endInd - k[2] <= nObs) {
      endInd <- min(endInd, nObs)

      # this part is only for approximating phi_t
      foo <- delta * allprobs[startInd, , drop = FALSE]
      foo <- foo / sum(foo)
      for(t in (startInd + 1):current_t) {
        foo <- (foo %*% Gamma) * allprobs[t, , drop = FALSE]
        foo <- foo / sum(foo)
      }
      # now we have phi_t and compute the likelihood contribution
      # from t+1 to t+k[2]
      for(t in (current_t+1):endInd) {
        foo <- (foo %*% Gamma) * allprobs[t, , drop = FALSE]
        sumfoo <- sum(foo)
        l <- l + log(sumfoo)
        foo <- foo / sumfoo
      }

      # updating indices
      current_t <- current_t + k[2]
      startInd <- current_t - k[1] + 1
      endInd <- current_t + k[2]
    }

  } else {
    delta1 <- delta1 + 1e-8
    delta1 <- delta1 / sum(delta1) # to avoid numerical issues

    delta <- delta + 1e-8
    delta <- delta / sum(delta) # to avoid numerical issues

    logdelta1 <- matrix(log(delta1), nrow = 1)
    logdelta <- matrix(log(delta), nrow = 1)

    logfoo <- logdelta1 + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:sum(k)){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l <- l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }

    current_t <- sum(k)
    startInd <- current_t - k[1] + 1
    endInd <- current_t + k[2]

    while(endInd - k[2] < nObs) {
      endInd <- min(endInd, nObs)
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      logfoo <- logdelta + logallprobs[startInd, , drop = FALSE]
      logfoo <- logfoo - LaMa:::logspace_add(logfoo)
      for(t in (startInd + 1):current_t) {
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        logfoo <- logfoo - LaMa:::logspace_add(logfoo)
      }
      for(t in (current_t+1):endInd) {
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        logsumfoo <- LaMa:::logspace_add(logfoo)
        l <- l + logsumfoo
        logfoo <- logfoo - logsumfoo
      }

      current_t <- current_t + k[2]
      startInd <- current_t - k[1] + 1
      endInd <- current_t + k[2]
    }
  }
  # TapeConfig(matmul = "atomic") # restore

  # returning approximate log-likelihood
  return(l)
}

# vector version
forward_banded3 <- function(delta, # initial distribution
                            Gamma, # tpm array
                            allprobs, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                            k = 20, # truncation parameter for sparse likelihood approximation
                            report = TRUE, # report quantities of interest?
                            logspace = FALSE # are allprobs on log-scale? -> internal computations in logspace
){
  # AD overloading (currently necessary with RTMB)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  # TapeConfig(matmul = "compact")

  if(logspace){
    logallprobs <- allprobs # assuming allprobs is on log-scale
    allprobs <- exp(logallprobs) # to report on the original scale
  }

  if(report){
    # reporting
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
  }

  if(length(k) == 1) k <- rep(k, 2)
  if(k[2] > k[1]) stop("k[2] must be smaller than or equal to k[1].")
  if(k[2] < 1) stop("k[2] must be at least 1.")

  # number of observations
  nObs <- nrow(allprobs)

  # assigning initial distribution as delta1
  delta1 <- delta
  # overwriting delta as stationary
  delta <- stationary(Gamma)

  # transposing Gamma
  Gammat = t(Gamma)

  if(!logspace){
    # initialise: compute f(x_1, ..., x_(k[1] + k[2]))
    # regular forward algorithm:
    foo <- delta1 * allprobs[1, ]
    sumfoo <- sum(foo)
    l <- log(sumfoo)
    foo <- foo / sumfoo
    for(t in 2:sum(k)){
      foo <- as.vector(Gammat %*% foo) * allprobs[t, ]
      sumfoo <- sum(foo)
      l <- l + log(sumfoo)
      foo <- foo / sumfoo
    }

    # now loop over tracks
    current_t <- sum(k) # t
    startInd <- current_t - k[1] + 1 # t - k[1] + 1
    endInd <- current_t + k[2] # t + k[2]

    while(current_t < nObs) {
      endInd <- min(endInd, nObs)

      # this part is only for approximating phi_t
      foo <- delta * allprobs[startInd, ]
      foo <- foo / sum(foo)
      for(t in (startInd + 1):current_t) {
        foo <- as.vector(Gammat %*% foo)  * allprobs[t, ]
        foo <- foo / sum(foo)
      }
      # now we have phi_t and compute the likelihood contribution
      # from t+1 to t+k[2]
      for(t in (current_t+1):endInd) {
        foo <- as.vector(Gammat %*% foo) * allprobs[t, ]
        sumfoo <- sum(foo)
        l <- l + log(sumfoo)
        foo <- foo / sumfoo
      }

      # updating indices
      current_t <- current_t + k[2]
      startInd <- current_t - k[1] + 1
      endInd <- current_t + k[2]
    }

  } else {
    delta1 <- delta1 + 1e-8
    delta1 <- delta1 / sum(delta1) # to avoid numerical issues

    delta <- delta + 1e-8
    delta <- delta / sum(delta) # to avoid numerical issues

    logdelta1 <- matrix(log(delta1), nrow = 1)
    logdelta <- matrix(log(delta), nrow = 1)

    logfoo <- logdelta1 + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:sum(k)){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l <- l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }

    current_t <- sum(k)
    startInd <- current_t - k[1] + 1
    endInd <- current_t + k[2]

    while(endInd - k[2] < nObs) {
      endInd <- min(endInd, nObs)
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      logfoo <- logdelta + logallprobs[startInd, , drop = FALSE]
      logfoo <- logfoo - LaMa:::logspace_add(logfoo)
      for(t in (startInd + 1):current_t) {
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        logfoo <- logfoo - LaMa:::logspace_add(logfoo)
      }
      for(t in (current_t+1):endInd) {
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        logsumfoo <- LaMa:::logspace_add(logfoo)
        l <- l + logsumfoo
        logfoo <- logfoo - logsumfoo
      }

      current_t <- current_t + k[2]
      startInd <- current_t - k[1] + 1
      endInd <- current_t + k[2]
    }
  }
  # TapeConfig(matmul = "atomic") # restore

  # returning approximate log-likelihood
  return(l)
}


# vector version
forward_banded4 <- function(delta, # initial distribution
                            Gamma, # tpm array
                            mu, sigma, gamma, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                            k = 20, # truncation parameter for sparse likelihood approximation
                            report = TRUE, # report quantities of interest?
                            logspace = FALSE # are allprobs on log-scale? -> internal computations in logspace
){
  # AD overloading (currently necessary with RTMB)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  if(report){
    # reporting
    REPORT(delta)
    REPORT(Gamma)
  }

  if(length(k) == 1) k <- rep(k, 2)
  if(k[2] > k[1]) stop("k[2] must be smaller than or equal to k[1].")
  if(k[2] < 1) stop("k[2] must be at least 1.")

  # number of observations
  nObs <- nrow(mu)

  # assigning initial distribution as delta1
  delta1 <- delta
  # overwriting delta as stationary
  delta <- stationary(Gamma)

  # transposing Gamma
  Gammat = t(Gamma)

  if(!logspace){
    # initialise: compute f(x_1, ..., x_(k[1] + k[2]))
    # regular forward algorithm:
    # calculate allprobs1
    allprobs <- dnorm(mu[1,1], 0, sigma) * dnorm(mu[1,2], 0, sigma)
    foo <- delta1 * allprobs
    sumfoo <- sum(foo)
    l <- log(sumfoo)
    foo <- foo / sumfoo
    for(t in 2:sum(k)){
      # calculate state-dep probs
      if(t == 2) {
        allprobs <- dnorm(mu[2,1], mu[1,1], sigma) * dnorm(mu[2,2], mu[1,2], sigma)
      } else {
        allprobs <- numeric(2)
        for(j in 1:2) {
          mean <- mu[t-1,] + gamma[j] * (mu[t-1,] - mu[t-2,])
          allprobs[j] <- dnorm(mu[t, 1], mean[1], sigma[j]) *
            dnorm(mu[t, 2], mean[2], sigma[j])
        }
      }

      foo <- as.vector(Gammat %*% foo) * allprobs
      sumfoo <- sum(foo)
      l <- l + log(sumfoo)
      foo <- foo / sumfoo
    }

    # now loop over tracks
    current_t <- sum(k) # t
    startInd <- current_t - k[1] + 1 # t - k[1] + 1
    endInd <- current_t + k[2] # t + k[2]

    while(current_t < nObs) {
      endInd <- min(endInd, nObs)

      allprobs <- numeric(2)
      for(j in 1:2) {
        mean <- mu[startInd-1,] + gamma[j] * (mu[startInd-1,] - mu[startInd-2,])
        allprobs[j] <- dnorm(mu[startInd, 1], mean[1], sigma[j]) *
          dnorm(mu[startInd, 2], mean[2], sigma[j])
      }

      # this part is only for approximating phi_t
      foo <- delta * allprobs
      foo <- foo / sum(foo)
      for(t in (startInd + 1):current_t) {
        allprobs <- numeric(2)
        for(j in 1:2) {
          mean <- mu[t-1,] + gamma[j] * (mu[t-1,] - mu[t-2,])
          allprobs[j] <- dnorm(mu[t, 1], mean[1], sigma[j]) *
            dnorm(mu[t, 2], mean[2], sigma[j])
        }

        foo <- as.vector(Gammat %*% foo)  * allprobs
        foo <- foo / sum(foo)
      }
      # now we have phi_t and compute the likelihood contribution
      # from t+1 to t+k[2]
      for(t in (current_t+1):endInd) {
        allprobs <- numeric(2)
        for(j in 1:2) {
          mean <- mu[t-1,] + gamma[j] * (mu[t-1,] - mu[t-2,])
          allprobs[j] <- dnorm(mu[t, 1], mean[1], sigma[j]) *
            dnorm(mu[t, 2], mean[2], sigma[j])
        }

        foo <- as.vector(Gammat %*% foo) * allprobs
        sumfoo <- sum(foo)
        l <- l + log(sumfoo)
        foo <- foo / sumfoo
      }

      # updating indices
      current_t <- current_t + k[2]
      startInd <- current_t - k[1] + 1
      endInd <- current_t + k[2]
    }

  } else {
    delta1 <- delta1 + 1e-8
    delta1 <- delta1 / sum(delta1) # to avoid numerical issues

    delta <- delta + 1e-8
    delta <- delta / sum(delta) # to avoid numerical issues

    logdelta1 <- matrix(log(delta1), nrow = 1)
    logdelta <- matrix(log(delta), nrow = 1)

    logfoo <- logdelta1 + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:sum(k)){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l <- l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }

    current_t <- sum(k)
    startInd <- current_t - k[1] + 1
    endInd <- current_t + k[2]

    while(endInd - k[2] < nObs) {
      endInd <- min(endInd, nObs)
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      logfoo <- logdelta + logallprobs[startInd, , drop = FALSE]
      logfoo <- logfoo - LaMa:::logspace_add(logfoo)
      for(t in (startInd + 1):current_t) {
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        logfoo <- logfoo - LaMa:::logspace_add(logfoo)
      }
      for(t in (current_t+1):endInd) {
        logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
        logsumfoo <- LaMa:::logspace_add(logfoo)
        l <- l + logsumfoo
        logfoo <- logfoo - logsumfoo
      }

      current_t <- current_t + k[2]
      startInd <- current_t - k[1] + 1
      endInd <- current_t + k[2]
    }
  }
  # TapeConfig(matmul = "atomic") # restore

  # returning approximate log-likelihood
  return(l)
}

# forward_track <- function(delta, Gamma, allprobs_track, k) {
#   end <- nrow(allprobs_track)
#   foo <- delta * allprobs_track[1, , drop = FALSE]
#   sumfoo <- sum(foo)
#   log_num <- log(sumfoo)
#   phi <- foo / sumfoo
#   for(t in 2:end){
#     foo <- (phi %*% Gamma) * allprobs_track[t, , drop = FALSE]
#     sumfoo <- sum(foo)
#     log_num <- log_num + log(sumfoo)
#     phi <- foo / sumfoo
#     # at the middle index, save value for denominator
#     if(t == k){
#       log_denom <- log_num
#     }
#   }
#   log_num - log_denom
# }
forward_track <- function(delta, Gamma, allprobs_track, k, logspace = FALSE) {
  end <- nrow(allprobs_track)
  track_l <- 0

  if(!logspace) {
    foo <- delta * allprobs_track[1, , drop = FALSE]
    foo <- foo / sum(foo)
    for(t in 2:k){
      foo <- (foo %*% Gamma) * allprobs_track[t, , drop = FALSE]
      foo <- foo / sum(foo)
    }
    for(t in (k+1):end){
      foo <- (foo %*% Gamma) * allprobs_track[t, , drop = FALSE]
      sumfoo <- sum(foo)
      track_l <- track_l + log(sumfoo)
      foo <- foo / sumfoo
    }
  } else{
    logdelta <- delta
    logallprobs_track <- allprobs_track

    logfoo <- logdelta + logallprobs_track[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    logfoo <- logfoo - logsumfoo
    for(t in 2:k){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs_track[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      logfoo <- logfoo - logsumfoo
    }
    for(t in (k+1):end){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs_track[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      track_l <- track_l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }
  }
  track_l
}

## banded version of the forward algorithm with one tpm
forward_banded2 <- function(delta, # initial distribution
                            Gamma, # tpm array
                            allprobs, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                            k = 30, # truncation parameter for sparse likelihood approximation
                            report = TRUE, # report quantities of interest?
                            logspace = FALSE # are allprobs on log-scale? -> internal computations in logspace
){
  # AD overloading (currently necessary with RTMB)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  if(logspace){
    logallprobs <- allprobs # assuming allprobs is on log-scale
    allprobs <- exp(logallprobs) # to report on the original scale
  }

  if(report){
    # reporting
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
  }

  # number of observations
  nObs <- nrow(allprobs)

  # assigning initial distribution as delta1
  # making delta1 a matrix for compatibility with RTMB
  delta1 <- matrix(delta, nrow = 1, ncol = length(delta))

  # overwriting delta as stationary
  delta <- stationary(Gamma)
  # making delta a matrix for compatibility with RTMB
  delta <- matrix(delta, nrow = 1, ncol = length(delta))

  n_tracks = ceiling(nObs / k) - 2 # number of remaining tracks
  trackInds <- list()
  startInd <- k + 1
  endInd <- k + 2*k
  for(track in seq_len(n_tracks)) {
    trackInds[[track]] <- startInd:(min(endInd, nObs))
    startInd <- startInd + k
    endInd <- endInd + k
  }

  if(!logspace){
    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    foo <- delta1 * allprobs[1, , drop = FALSE]
    sumfoo <- sum(foo)
    l <- log(sumfoo)
    phi <- foo / sumfoo
    for(t in 2:(2*k)){
      foo <- (phi %*% Gamma) * allprobs[t, , drop = FALSE]
      # foo <- (phi %sp% Gamma) * allprobs[t, , drop = FALSE]
      sumfoo <- sum(foo)
      l <- l + log(sumfoo)
      phi <- foo / sumfoo
    }

    # now we always compute f(x_t, ... x_t+k | x_t-1, ..., x_t-k)
    # = f(x_t-k, ..., x_t-1, x_t, ... x_t+k) / f(x_t-k, ..., x_t-1)
    # the nice thing is that both the numerator and denominator can be
    # computed with one pass of the forward algorithm
    lls <- sapply(1:n_tracks, function(track){
      forward_track(delta, Gamma, allprobs[trackInds[[track]], , drop = FALSE], k)
    })
    return(l + sum(lls))
  } else {
    delta1 <- delta1 + 1e-8
    delta1 <- delta1 / sum(delta1) # to avoid numerical issues

    delta <- delta + 1e-8
    delta <- delta / sum(delta) # to avoid numerical issues

    logdelta1 <- matrix(log(delta1), nrow = 1)
    logdelta <- matrix(log(delta), nrow = 1)

    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    logfoo <- logdelta1 + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:(2*k)){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
      # logfoo <- log(exp(logfoo) %sp% Gamma) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l <- l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }

    # now we always compute f(x_t, ... x_t+k | x_t-1, ..., x_t-k)
    # = f(x_t-k, ..., x_t-1, x_t, ... x_t+k) / f(x_t-k, ..., x_t-1)
    # the nice thing is that both the numerator and denominator can be
    # computed with one pass of the forward algorithm
    lls <- sapply(1:n_tracks, function(track){
      forward_track(logdelta, Gamma, logallprobs[trackInds[[track]], , drop = FALSE], k, logspace = TRUE)
    })
    return(l + sum(lls))
  }
}



## banded version of the forward algorithm with tpm array
forward_g_banded <- function(delta, # initial distribution
                             Gamma, # tpm array
                             allprobs, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                             k = 30, # truncation parameter for sparse likelihood approximation
                             report = TRUE, # report quantities of interest?
                             logspace = TRUE, # are allprobs on log-scale? -> internal computations in logspace
                             j = k # truncation parameter for calculation of state distribution, if NULL, no approximation is used
                             ){

  # AD overloading (currently necessary with RTMB)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  if(logspace){
    logallprobs <- allprobs # assuming allprobs is on log-scale
    allprobs <- exp(logallprobs) # to report on the original scale
  }

  if(report){
    # reporting
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
  }

  # number of observations
  nObs <- nrow(allprobs)

  if (dim(Gamma)[3] == nObs - 1) {
    Gamma <- AD(array(c(array(0, dim = dim(Gamma)[1:2]), Gamma),
                      dim = c(dim(Gamma)[1], dim(Gamma)[2], nObs)))
  }

  # number of remaining tracks after first two exact tracks
  n_tracks = ceiling(nObs / k) - 2

  # where does the last track start? We only need state distribution until here
  last_start = n_tracks * k + k + 1

  # computing all unconditional state distributions at each time point
  if(is.null(j)){
    Delta <- statedist(delta, Gamma[, , seq_len(last_start)])
  } else {
    Delta <- statedist_banded(delta, Gamma[, , seq_len(last_start)], j = j)
  }

  # initialise log-likelihood value
  l <- 0

  if(!logspace){
    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    foo <- Delta[1, , drop = FALSE] * allprobs[1, , drop = FALSE]
    sumfoo <- sum(foo)
    l <- log(sumfoo)
    phi <- foo / sumfoo
    for(t in 2:(2*k)){
      foo <- (phi %*% Gamma[,, t]) * allprobs[t, , drop = FALSE]
      sumfoo <- sum(foo)
      l <- l + log(sumfoo)
      phi <- foo / sumfoo
    }

    # now we always compute f(x_t, ... x_t+k | x_t-1, ..., x_t-k)
    # = f(x_t-k, ..., x_t-1, x_t, ... x_t+k) / f(x_t-k, ..., x_t-1)
    # the nice thing is that both the numerator and denominator can be
    # computed with one pass of the forward algorithm
    startInd <- k + 1
    endInd <- k + 2*k
    for(track in seq_len(n_tracks)){
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      foo <- Delta[startInd, , drop = FALSE] * allprobs[startInd, , drop = FALSE]
      sumfoo <- sum(foo)
      log_num <- log(sumfoo)
      phi <- foo / sumfoo
      end <- min(endInd, nObs)
      for(t in (startInd + 1):end){
        foo <- (phi %*% Gamma[,, t]) * allprobs[t, , drop = FALSE]
        sumfoo <- sum(foo)
        log_num <- log_num + log(sumfoo)
        phi <- foo / sumfoo
        # at the middle index, save value for denominator
        if(t == (startInd + k - 1)){
          log_denom <- log_num
        }
      }
      l <- l + log_num - log_denom

      startInd <- startInd + k
      endInd <- endInd + k
    }
  } else {
    Delta[1, ] <- Delta[1,] + 1e-8
    Delta[1, ] <- Delta[1, ] / sum(Delta[1, ]) # to avoid numerical issues w log(0)

    logDelta <- log(Delta)

    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    logfoo <- logDelta[1, , drop = FALSE] + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:(2*k)){
      logfoo <- log(exp(logfoo) %*% Gamma[,, t]) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l <- l + logsumfoo
      logfoo <- logfoo - logsumfoo
    }

    # now we always compute f(x_t, ... x_t+k | x_t-1, ..., x_t-k)
    # = f(x_t-k, ..., x_t-1, x_t, ... x_t+k) / f(x_t-k, ..., x_t-1)
    # the nice thing is that both the numerator and denominator can be
    # computed with one pass of the forward algorithm
    startInd <- k + 1
    endInd <- k + 2*k
    for(track in seq_len(n_tracks)){
      # regular forward algorithm
      # what creates sparsity/ banded hessian here is of course that we start with delta,
      # not with phi for the conditioning
      logfoo <- logDelta[startInd, , drop = FALSE] + logallprobs[startInd, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      log_num <- logsumfoo
      logfoo <- logfoo - logsumfoo
      end <- min(endInd, nObs)
      for(t in (startInd + 1):end){
        logfoo <- log(exp(logfoo) %*% Gamma[,, t]) + logallprobs[t, , drop = FALSE]
        logsumfoo <- LaMa:::logspace_add(logfoo)
        log_num <- log_num + logsumfoo
        logfoo <- logfoo - logsumfoo
        # at the middle index, save value for denominator
        if(t == (startInd + k - 1)){
          log_denom <- log_num
        }
      }
      l <- l + log_num - log_denom

      startInd <- startInd + k
      endInd <- endInd + k
    }
  }

  # returning approximate log-likelihood
  l
}








# Compute k ---------------------------------------------------------------

forward_summands <- function(delta, # initial distribution
                           Gamma, # tpm array
                           allprobs, # matrix of state-dependent densities (on logscale if logspace = TRUE)
                           logspace = FALSE # are allprobs on log-scale? -> internal computations in logspace
){

  if(logspace){
    logallprobs <- allprobs # assuming allprobs is on log-scale
  }

  # number of observations
  nObs <- nrow(allprobs)

  # assigning initial distribution as delta1
  # making delta1 a matrix for compatibility with RTMB
  delta1 <- matrix(delta, nrow = 1, ncol = length(delta))

  # overwriting delta as stationary
  delta <- stationary(Gamma)
  # making delta a matrix for compatibility with RTMB
  delta <- matrix(delta, nrow = 1, ncol = length(delta))

  l <- rep(NA, nObs)

  if(!logspace){
    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    foo <- delta1 * allprobs[1, , drop = FALSE]
    sumfoo <- sum(foo)
    l[1] <- log(sumfoo)
    phi <- foo / sumfoo
    for(t in 2:nObs){
      foo <- (phi %*% Gamma) * allprobs[t, , drop = FALSE]
      # foo <- (phi %sp% Gamma) * allprobs[t, , drop = FALSE]
      sumfoo <- sum(foo)
      l[t] <- log(sumfoo)
      phi <- foo / sumfoo
    }

  } else {
    delta1 <- delta1 + 1e-8
    delta1 <- delta1 / sum(delta1) # to avoid numerical issues

    delta <- delta + 1e-8
    delta <- delta / sum(delta) # to avoid numerical issues

    logdelta1 <- matrix(log(delta1), nrow = 1)
    logdelta <- matrix(log(delta), nrow = 1)

    # initialise: compute f(x_1, ..., x_2k)
    # regular forward algorithm:
    logfoo <- logdelta1 + logallprobs[1, , drop = FALSE]
    logsumfoo <- LaMa:::logspace_add(logfoo)
    l[1] <- logsumfoo
    logfoo <- logfoo - logsumfoo
    for(t in 2:nObs){
      logfoo <- log(exp(logfoo) %*% Gamma) + logallprobs[t, , drop = FALSE]
      # logfoo <- log(exp(logfoo) %sp% Gamma) + logallprobs[t, , drop = FALSE]
      logsumfoo <- LaMa:::logspace_add(logfoo)
      l[t] <- logsumfoo
      logfoo <- logfoo - logsumfoo
    }
  }

  return(l)
}
