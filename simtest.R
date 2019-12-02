library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")

# test case 1
T <- 100
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.2,0.75),K,K,byrow = TRUE)
B <- matrix(c(0.5,0.1,0.3,0.4,0.4,0.4,0.2,0.6,0.2),K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

A <- list()
for(i in 1:m){
  edge <- sample(m, 25) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))


# just test how fast simulation is 
#start_time <- Sys.time()
#for(i in 1:100000){
#  x <- sampleChild(0,100,1,0.5)
#}
#end_time <- Sys.time()
#print(end_time - start_time)

# main
dT <- 0.5
tau <- matrix(0,m,K)
for (k in 1:K){
  tau[which(Z == (k-1)),k] <- 1
}
system.time(results.online <- online_estimator(alltimes, A, m, K, T, dT, 
                                               lam = 1.0, B, Mu, tau))

system.time(results.eff <- online_estimator_eff(alltimes, A, m, K, T, dT,
                                                lam = 1.0, B, Mu, tau)) 

H <- 1
window <- 1
MuA <- array(0,c(K,K,H))
system.time(results.nh.eff <- nonhomoHak_estimator_eff(alltimes,A,m,K,H,
                                                    window,T,dT,lam = 0.1, gravity = 0.0, B,MuA,tau))


#system.time(results.stoch <- online_estimator(alltimes, A, m, K, T, dT, 
#                                               lam = 1.0, B, Mu, tau, percent = 0.25))

itermax <- T / dT
stop_eps <- 0.001
system.time(results.batch <- batch_estimator(alltimes, A, m, K, T, dT, 
                                             lam = 1.0, B, Mu, tau, itermax, stop_eps))


# -- calculate loglik ---
Z <- apply(results.online$tau,1,which.max) - 1
loglik.online <- get_loglik_Hak(alltimes,0,T,Z,results.online$Mu, 
                                results.online$B, results.online$Pi,A, results.online$lam,m,K)

Z <- apply(results.eff$tau,1,which.max) - 1
loglik.eff <- get_loglik_Hak(alltimes,0,T,Z,results.eff$Mu, 
                                results.eff$B, results.eff$Pi,A, results.eff$lam,m,K)


Z <- apply(results.batch$tau,1,which.max) - 1
loglik.batch <- get_loglik_Hak(alltimes,0,T,Z,results.batch$Mu, 
                                results.batch$B, results.batch$Pi,A, results.batch$lam,m,K)


1 - (loglik.batch - loglik.eff)/abs(loglik.batch)



# -- debug elbo ---


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
