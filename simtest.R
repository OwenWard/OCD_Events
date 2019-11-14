# test case 1
T <- 100
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


# -- debug elbo ---
results1 <- results
results2 <- results



tau <- results1$tau
Mu <- results1$B
Pi <- results1$Pi
lam <- results1$lam
elbo1 <- get_elbo_hak(alltimes,0,T,tau,Mu,B,Pi,A,lam,m,K)

results2 <- results
tau <- results2$tau
Mu <- results2$B
Pi <- results2$Pi
lam <- results2$lam
elbo2 <- get_elbo_hak(alltimes,0,T,tau,Mu,B,Pi,A,lam = 1.0,m,K)

# -- debug lam update---
# grad_lam <- test_lam(tau, Mu, B, Pi, S, alltimes, 0, T, m, K, A, lam = 0.9)



### Poisson Simulation ###
T = 100
dT = 0.5
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B <- matrix(0,K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))

# initial
Pi = c(0.3,0.3,0.4)
B = matrix(c(1.2,0.5,0.5,0.5,1.1,.65,0.75,0.85,1.15),nrow = K,ncol = K,byrow = T)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results <- estimate_Poisson(full_data = alltimes,tau,B,Pi,S,A,m,K,dT,T)
results$B
