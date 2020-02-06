library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")

# test case 1
T <- 200
K <- 1
b <- 0.5
m <- 5
W1 <- matrix(0,m,K)
W1[,1] <- c(rep(1,m))
W2 <- matrix(0,m,K)
W2[,1] <- c(rep(1,m))

A <- list()
for(i in 1:m){
  edge <- sample(m, 5) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


system.time(alltimes <- sampleCCRMHak(T,A,W1,W2,b,lam = 1))

dT <- 1.0
system.time(results.ccrm <- ccrm_estimator(alltimes,A,m,K,T,dT,lam = 1.0,W1,W2,b))
system.time(results.ccrmbatch <- batch_ccrm_estimator(alltimes,A,m,K,T,dT,lam = 1.0,
                                                      W1,W2,b,50,0.001))





dT <- 1
W <- matrix(0,m,m)
system.time(results.test <- test_estimator(alltimes,A,m,K,T,dT,lam = 1.0,W,b))


# -----------------------
# test case 2
T <- 200
K <- 2
b <- 0.5
m <- 10
W1 <- matrix(0,m,K)
W1[,1] <- c(rep(1,m/2),rep(0.1,m/2))
W1[,2] <- c(rep(0.1,m/2),rep(1,m/2))
W2 <- matrix(0,m,K)
W2[,1] <- c(rep(1,m/2),rep(0.1,m/2))
W2[,2] <- c(rep(0.1,m/2),rep(1,m/2))

A <- list()
for(i in 1:m){
  edge <- sample(m, m) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


system.time(alltimes <- sampleCCRMHak(T,A,W1,W2,b,lam = 1))

dT <- 1.0
system.time(results.ccrm <- ccrm_estimator(alltimes,A,m,K,T,dT,lam = 1.0,W1,W2,b))



# -----------------------
# test case 2
T <- 200
K <- 3
b <- 0.5
m <- 50
W1 <- matrix(0,m,K)
W1[,1] <- c(rep(0.8,m/2),rep(0.1,m/2))
W1[,2] <- c(rep(0.1,m/2),rep(0.8,m/2))
W1[,3] <- c(rep(0.4,m/2),rep(0.4,m/2))
W2 <- matrix(0,m,K)
W2[,1] <- c(rep(0.8,m/2),rep(0.1,m/2))
W2[,2] <- c(rep(0.1,m/2),rep(0.8,m/2))
W2[,3] <- c(rep(0.4,m/2),rep(0.4,m/2))

A <- list()
for(i in 1:m){
  edge <- sample(m, 15) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


system.time(alltimes <- sampleCCRMHak(T,A,W1,W2,b,lam = 1))

dT <- 1.0
system.time(results.ccrm <- ccrm_estimator(alltimes,A,m,K,T,dT,lam = 1.0,W1,W2,b))


