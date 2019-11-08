# test case 1
T <- 50
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B <- matrix(c(0.5,0.1,0.3,0.4,0.4,0.4,0.2,0.3,0.6),K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

A <- list()
for(i in 1:m){
  edge <- sample(m, 25) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))

# main
dT <- 0.25
tau <- matrix(0,m,K)
for (k in 1:K){
  tau[which(Z == (k-1)),k] <- 1
}
system.time(results <- online_estimator(alltimes, A, m, K, T, dT, lam = 1, B, Mu, tau))


# --- test ----
S <- matrix(0,m,K)
tau <- matrix(0,m,K)
for (k in 1:K){
  tau[which(Z == (k-1)),k] <- 1
}
paralist <- test(alltimes,tau,Mu,B,Pi,S,0,50,m,K,A,1,1/dim(alltimes)[1])


