library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")

# test case 1
T <- 100
K <- 2
H <- 2
MuA <- array(0, c(K, K, H))
MuA[, , 1] <- matrix(c(0.8, 0.2, 0.6, 0.4), 2, 2)
MuA[, , 2] <- matrix(c(0.4, 0.7, 0.2, 0.7), 2, 2)
B <- matrix(c(0.8, 0.2, 0.4, 0.6), K, K, byrow = TRUE)
m <- 100
Pi <- matrix(c(0.6, 0.4), 1, K)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]))
window <- 0.25

A <- list()
for (i in 1:m) {
  edge <- sample(m, 25) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}



system.time(alltimes <- sampleBlockHak_nonhomo(T, A, Z, MuA, B, window, lam = 1))

dT <- 0.5
tau <- matrix(0, m, K)
for (k in 1:K) {
  tau[which(Z == (k - 1)), k] <- 1
}

system.time(results.online <- nonhomoHak_estimator(alltimes,
  A,
  m,
  K,
  H,
  window,
  T,
  dT,
  lam = 0.1,
  gravity = 0.0,
  B,
  MuA,
  tau
))

system.time(results.eff <- nonhomoHak_estimator_eff(alltimes, A, m, K, H,
  window, T, dT,
  lam = 0.1, gravity = 0.0, B, MuA, tau
))

itermax <- T / dT
stop_eps <- 0.001
system.time(results.batch <- batch_nonhomoHak_estimator(alltimes, A, m, K, H,
  window, T, dT,
  lam = 0.1, gravity = 0.0,
  B, MuA, tau, itermax, stop_eps
))

# -- calculate loglik ---
Z <- apply(results.eff$tau, 1, which.max) - 1
loglik.eff <- get_loglik_nonhomoHak(
  alltimes, 0, T, Z, results.eff$MuA,
  results.eff$B, results.eff$Pi, A, results.eff$lam,
  m, K, H, window
)

Z <- apply(results.batch$tau, 1, which.max) - 1
loglik.batch <- get_loglik_nonhomoHak(
  alltimes, 0, T, Z, results.batch$MuA,
  results.batch$B, results.batch$Pi, A, results.batch$lam,
  m, K, H, window
)

1 - (loglik.batch - loglik.eff) / abs(loglik.batch)

# -- debug elbo ---
# elbo <- get_elbo_nonhomoHak(alltimes,0,T,tau,MuA,B,Pi,A,lam = 0.5,m,K,H,window)




# test case 2
T <- 500
K <- 2
H <- 7
MuA <- array(0, c(K, K, H))
MuA[, , 1] <- matrix(c(0.8, 0.2, 0.3, 0.6), 2, 2, byrow = T)
MuA[, , 2] <- matrix(c(0, 0.5, 0.2, 0.7), 2, 2, byrow = T)
MuA[, , 3] <- matrix(c(1.0, 0, 0.5, 0.6), 2, 2, byrow = T)
MuA[, , 4] <- matrix(c(0.4, 0.4, 0, 0.7), 2, 2, byrow = T)
MuA[, , 5] <- matrix(c(0.5, 0.5, 0.5, 0), 2, 2, byrow = T)
MuA[, , 6] <- matrix(c(0.4, 0, 0.4, 0), 2, 2, byrow = T)
MuA[, , 7] <- matrix(c(0, 0.9, 0, 0.9), 2, 2, byrow = T)
B <- matrix(c(0.8, 0.3, 0.4, 0.7), K, K, byrow = TRUE)
m <- 500
Pi <- matrix(c(0.6, 0.4), 1, K)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]))
window <- 1 / 7

A <- list()
for (i in 1:m) {
  edge <- sample(m, 25) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak_nonhomo(T, A, Z, MuA, B, window, lam = 1))

dT <- 2.5
tau <- matrix(0, m, K)
for (k in 1:K) {
  tau[which(Z == (k - 1)), k] <- 1
}

system.time(results.online <- nonhomoHak_estimator(alltimes, A, m, K, H, window, T, dT,
  lam = 0.1,
  gravity = 0.001, B, MuA, tau
))

itermax <- T / dT
stop_eps <- 0.001
system.time(results.batch <- batch_nonhomoHak_estimator(alltimes, A, m, K, H,
  window, T, dT,
  lam = 0.1,
  gravity = 0.001,
  B, MuA, tau, itermax, stop_eps
))

# -- calculate loglik ---
Z <- apply(results.online$tau, 1, which.max) - 1
loglik.online <- get_loglik_nonhomoHak(
  alltimes, 0, T, Z, results.online$MuA,
  results.online$B, results.online$Pi, A, results.online$lam,
  m, K, H, window
)

Z <- apply(results.batch$tau, 1, which.max) - 1
loglik.batch <- get_loglik_nonhomoHak(
  alltimes, 0, T, Z, results.batch$MuA,
  results.batch$B, results.batch$Pi, A, results.batch$lam,
  m, K, H, window
)

1 - (loglik.batch - loglik.online) / abs(loglik.batch)


# test case 3 --- nonhomo poisson case
T <- 100
K <- 2
H <- 7
MuA <- array(0, c(K, K, H))
MuA[, , 1] <- matrix(c(0.8, 0.2, 0.3, 0.6), 2, 2, byrow = T)
MuA[, , 2] <- matrix(c(0, 0.5, 0.2, 0.7), 2, 2, byrow = T)
MuA[, , 3] <- matrix(c(1.0, 0, 0.5, 0.6), 2, 2, byrow = T)
MuA[, , 4] <- matrix(c(0.4, 0.4, 0, 0.7), 2, 2, byrow = T)
MuA[, , 5] <- matrix(c(0.5, 0.5, 0.5, 0), 2, 2, byrow = T)
MuA[, , 6] <- matrix(c(0.4, 0, 0.4, 0), 2, 2, byrow = T)
MuA[, , 7] <- matrix(c(0, 0.9, 0, 0.9), 2, 2, byrow = T)
B <- matrix(0, K, K, byrow = TRUE)
m <- 100
Pi <- matrix(c(0.6, 0.4), 1, K)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]))
window <- 1 / 7

A <- list()
for (i in 1:m) {
  edge <- sample(m, 25) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak_nonhomo(T, A, Z, MuA, B, window, lam = 1))

dT <- 0.5
tau <- matrix(0, m, K)
for (k in 1:K) {
  tau[which(Z == (k - 1)), k] <- 1
}

system.time(results.online <- nonhomoPois_estimator(alltimes, A, m, K, H, window, T, dT,
  gravity = 0.001, MuA, tau
))

itermax <- T / dT
stop_eps <- 0.001
system.time(results.batch <- batch_nonhomoPois_estimator(alltimes, A, m, K, H, window, T, dT,
  gravity = 0.001, MuA, tau, itermax, stop_eps
))


# -- calculate loglik ---
Z <- apply(results.online$tau, 1, which.max) - 1
loglik.online <- get_loglik_nonhomoHak(
  alltimes, 0, T, Z, results.online$MuA,
  B, results.online$Pi, A, 0.0,
  m, K, H, window
)

Z <- apply(results.batch$tau, 1, which.max) - 1
loglik.batch <- get_loglik_nonhomoHak(
  alltimes, 0, T, Z, results.batch$MuA,
  B, results.batch$Pi, A, 0.0,
  m, K, H, window
)

1 - (loglik.batch - loglik.online) / abs(loglik.batch)
