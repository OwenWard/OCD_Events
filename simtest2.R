# test case 1
T <- 50
K <- 2
H <- 2
MuA <- array(0,c(K,K,H))
MuA[,,1] <- matrix(c(0.8,0.2,0.6,0.4),2,2)
MuA[,,2] <- matrix(c(0.4,0.7,0.2,0.7),2,2)
B <- matrix(c(0.8,0.2,0.4,0.6),K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.6,0.4),1,K)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]))
window <- 0.25

A <- list()
for(i in 1:m){
  edge <- sample(m, 25) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")

system.time(alltimes <- sampleBlockHak_nonhomo(T, A, Z, MuA, B, window, lam = 1))

dT <- 0.25
tau <- matrix(0,m,K)
for (k in 1:K){
  tau[which(Z == (k-1)),k] <- 1
}

system.time(results <- nonhomoHak_estimator(alltimes,A,m,K,H,window,T,dT,lam = 0.1,B,MuA,tau))

# -- debug elbo ---
elbo <- get_elbo_nonhomoHak(alltimes,0,T,tau,MuA,B,Pi,A,lam = 0.5,m,K,H,window)

