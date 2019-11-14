# because these models don't scale, use on a subset of the data?
library(mclust)
library(dplyr)
library(ppsbm)

library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")


#### Simulate data from our Model, fit with all methods ####
T = 50
dT = 0.1
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B <- matrix(0,K,K,byrow = TRUE)
m <- 30
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

A <- list()
for(i in 1:m){
  edge <- sample(m, 10) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))

Pi = c(0.3,0.3,0.4)
B = matrix(c(1,0.5,0.5,0.5,1,.65,0.75,0.85,1.),nrow = K,ncol = K,byrow = T)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

# more random start
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results_online <- estimate_Poisson(full_data = alltimes,tau,B,Pi,S,A,m,K,dT,T)
results_online$B
Mu
est_Z = apply(results_online$tau,1,which.max)#-1



adjustedRandIndex(Z,est_Z)

# then transform and fit with PPSBM
alltimes %>% head()

n_events = dim(alltimes)[1]
ids = rep(0,n_events)
times = rep(0,n_events)

A_matrix = matrix(0,nrow = m,ncol = m)
A_mat_bin = matrix(0,nrow = m,ncol = m)

for(i in 1:n_events){
  j = alltimes[i,1]+1
  k = alltimes[i,2]+1
  ids[i] = convertNodePair(j,k,m,TRUE)
  times[i] = alltimes[3]
  A_matrix[j,k] = A_matrix[j,k] + 1
  A_mat_bin[j,k] = 1
}

sim_data = list(time.seq = times, type.seq = ids, Time = T)

results = mainVEM(data = sim_data, n = m, Qmin = 3, Qmax = 3, method = 'kernel',
                  d_part = 1, n_perturb = 1)

results[[2]]$tau
est_Z_ppsbm = apply(results[[1]]$tau, 2, which.max)
#est_Z_ppsbm
adjustedRandIndex(Z,est_Z_ppsbm)



### Spectral Clustering on Weighted Adjacency Matrix ####
A_matrix

library(kernlab)

sc <- specc(A_matrix, centers=3)
estimated_groups = sc@.Data
adjustedRandIndex(estimated_groups,Z)


### SBM on Binary Counts ####
library(blockmodels)
sbm_bin = BM_bernoulli('SBM',A_mat_bin,explore_min = 3, explore_max = 3)
sbm_bin$estimate()
sbm_bin$memberships

sbm_bin_est = apply(sbm_bin$memberships[[3]]$Z,1,which.max)
adjustedRandIndex(Z,sbm_bin_est)


### SBM on Counts ####
sbm_count = BM_poisson('SBM',A_matrix, explore_min = 3, explore_max = 3)
sbm_count$estimate()
sbm_count_est = apply(sbm_count$memberships[[3]]$Z,1,which.max)
adjustedRandIndex(Z,sbm_count_est)
