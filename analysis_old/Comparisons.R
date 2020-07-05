# because these models don't scale, use on a subset of the data?
# library(mclust)
# library(dplyr)
# library(ppsbm)
# library(kernlab)
# library(blockmodels)
# library(Rcpp)
# library(RcppArmadillo)
# sourceCpp("onlineblock.cpp")
source("comparison_fcns.R")



#### Simulate data from our Model, fit with all methods ####
T = 50
dT = 0.1
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B_Pois <- matrix(0,K,K,byrow = TRUE)
m <- 300
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  edge <- sample(m, 100) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B_Pois, lam = 1))

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
#alltimes %>% head()

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
                  directed = T,
                  d_part = 1, n_perturb = 1)

#results[[2]]$tau
est_Z_ppsbm = apply(results[[1]]$tau, 2, which.max)
#est_Z_ppsbm
adjustedRandIndex(Z,est_Z_ppsbm)

steps = 10
for(i in 1:10){
  est_Z = apply(results_online$early_tau[,,i],1,which.max)
  print(adjustedRandIndex(Z,est_Z))
}



### Spectral Clustering on Weighted Adjacency Matrix ####

spec_ARI(A_matrix,K=3,Z)

### SBM on Binary Counts ####


sbm_bin_ARI(A_mat_bin,K=3,Z)


### SBM on Counts ####

sbm_count_ARI(A_matrix,K=3,Z)



#### Simulate Data from Hawkes Process and Fit the Models again
Mu_H <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B_H <- matrix(c(0.7,0.1,0.1,0.15,0.85,0.2,0.01,0.4,1),K,K,byrow = TRUE)


system.time(alltimes_Hawkes <- sampleBlockHak(T, A, Z, Mu_H, B_H, lam = 1))
n_events = dim(alltimes_Hawkes)[1]
ids = rep(0,n_events)
times = rep(0,n_events)

A_matrix_H = matrix(0,nrow = m,ncol = m)
A_mat_bin_H = matrix(0,nrow = m,ncol = m)

for(i in 1:n_events){
  j = alltimes_Hawkes[i,1]+1
  k = alltimes_Hawkes[i,2]+1
  ids[i] = convertNodePair(j,k,m,TRUE)
  times[i] = alltimes[3]
  A_matrix_H[j,k] = A_matrix[j,k] + 1
  A_mat_bin_H[j,k] = 1
}

sim_data_H = list(time.seq = times, type.seq = ids, Time = T)

results_H = mainVEM(data = sim_data_H, n = m, Qmin = 3, Qmax = 3, method = 'kernel',
                  directed = T,
                  d_part = 1, n_perturb = 1)

est_Z_ppsbm_H = apply(results_H[[1]]$tau, 2, which.max)
#est_Z_ppsbm
adjustedRandIndex(Z,est_Z_ppsbm_H)


# fit our hawkes model also
Pi = c(0.3,0.3,0.4)
B = matrix(c(1,0.5,0.5,0.5,1,.65,0.75,0.85,1.),nrow = K,ncol = K,byrow = T)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)
Mu = matrix(runif(K*K),K,K)


system.time(results_H <- online_estimator(alltimes, A, m, K, T, dT, lam = 1, B, Mu, tau))
# true mu and B on lines 16 and 32
est_Z_H <- apply(results_H$tau, 1, which.max)
adjustedRandIndex(Z,est_Z_H)

# spectral clustering on Hawkes ####

spec_ARI(A_matrix_H,K=3,Z)

# Binary SBM #### 

sbm_bin_ARI(A_mat_bin_H,3,Z)

# Poisson SBM ###

sbm_count_ARI(A_matrix_H,3,Z)

##### Run a bunch of simulations Poisson data ####
n_iters = 50

Pois_On_ARI <- rep(0,n_iters)
Pois_PPSBM_ARI <- rep(0,n_iters)
Pois_Sc_ARI <- rep(0,n_iters)
Pois_B_SBM_ARI <- rep(0,n_iters)
Pois_P_SBM_ARI <- rep(0,n_iters)

for(l in 1:n_iters){
  
  # simulate data ####
  T = 50
  dT = 0.1
  K <- 3
  Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
  B_Pois <- matrix(0,K,K,byrow = TRUE)
  m <- 30
  Pi <- matrix(c(0.4,0.3,0.3),1,3)
  Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))
  
  A <- list()
  for(i in 1:m){
    edge <- sample(m, 10) - 1
    edge <- sort(edge)
    A[[i]] <- edge
  }
  
  system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B_Pois, lam = 1))
  
  ### fit our model ####
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
  Pois_On_ARI[l] = adjustedRandIndex(Z,est_Z)
  ### fit ppsbm ####
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
  
  # results = mainVEM(data = sim_data, n = m, Qmin = 3, Qmax = 3, method = 'kernel',
  #                   directed = T,
  #                   d_part = 1, n_perturb = 1)
  # 
  # #results[[2]]$tau
  # est_Z_ppsbm = apply(results[[1]]$tau, 2, which.max)
  # #est_Z_ppsbm
  # Pois_PPSBM_ARI[l] = adjustedRandIndex(Z,est_Z_ppsbm)
  
  Pois_Sc_ARI[l] = spec_ARI(A_matrix,K=3,Z)
  
  ### SBM on Binary Counts ####
  
  
  Pois_B_SBM_ARI[l] = sbm_bin_ARI(A_mat_bin,K=3,Z)
  
  
  ### SBM on Counts ####
  
  Pois_P_SBM_ARI[l] = sbm_count_ARI(A_matrix,K=3,Z)
  
  
}

# then compare these
mean(Pois_On_ARI)
mean(Pois_Sc_ARI)
mean(Pois_P_SBM_ARI)
mean(Pois_B_SBM_ARI)
# performance similar to Poisson SBM when Poisson data. our model doesn't use multiple 
# restarts, etc


### Simulate from Hawkes model ####
n_iters = 20

Hawkes_On_ARI <- rep(0,n_iters)
Pois_PPSBM_ARI <- rep(0,n_iters)
Pois_Sc_ARI <- rep(0,n_iters)
Pois_B_SBM_ARI <- rep(0,n_iters)
Pois_P_SBM_ARI <- rep(0,n_iters)

for(l in 1:n_iters){
  
  # simulate data ####
  T = 50
  dT = 0.2
  K <- 3
  Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
  B = matrix(c(1,0.5,0.5,0.5,1,.65,0.75,0.85,1.),nrow = K,ncol = K,byrow = T)
  m <- 50
  Pi <- matrix(c(0.4,0.3,0.3),1,3)
  Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))
  
  A <- list()
  for(i in 1:m){
    edge <- sample(m, 10) - 1
    edge <- sort(edge)
    A[[i]] <- edge
  }
  
  system.time(alltimes <- sampleBlockHak(T, A, Z, Mu, B, lam = 1))
  
  ### fit our model ####
  Pi = rep(1/K,K)
  B = matrix(runif(K*K),K,K)
  Mu = matrix(runif(K*K),K,K)
  diag(Mu) = rnorm(3,mean = 1, sd = 0.1)
  tau = matrix(runif(m*K),nrow=m,ncol=K)
  tau = tau/rowSums(tau)
  S = matrix(0,nrow = m,ncol = K)
  
  
  # results_online <- estimate_Poisson(full_data = alltimes,tau,B,Pi,S,A,m,K,dT,T)
  # results_online$B
  # Mu
  # est_Z = apply(results_online$tau,1,which.max)#-1
  # Pois_On_ARI[l] = adjustedRandIndex(Z,est_Z)
  system.time(results_H <- online_estimator(alltimes, A, m, K, T, dT, lam = 1, B, Mu, tau))
  # true mu and B on lines 16 and 32
  est_Z_H <- apply(results_H$tau, 1, which.max)
  Hawkes_On_ARI[l] = adjustedRandIndex(Z,est_Z_H)
  
  
  
  ### fit ppsbm ####
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
  
  # results = mainVEM(data = sim_data, n = m, Qmin = 3, Qmax = 3, method = 'kernel',
  #                   directed = T,
  #                   d_part = 1, n_perturb = 1)
  # 
  # #results[[2]]$tau
  # est_Z_ppsbm = apply(results[[1]]$tau, 2, which.max)
  # #est_Z_ppsbm
  # Pois_PPSBM_ARI[l] = adjustedRandIndex(Z,est_Z_ppsbm)
  
  Pois_Sc_ARI[l] = spec_ARI(A_matrix,K=3,Z)
  
  ### SBM on Binary Counts ####
  
  
  Pois_B_SBM_ARI[l] = sbm_bin_ARI(A_mat_bin,K=3,Z)
  
  
  ### SBM on Counts ####
  
  Pois_P_SBM_ARI[l] = sbm_count_ARI(A_matrix,K=3,Z)
  
  
}

# then compare these
mean(Hawkes_On_ARI)
mean(Pois_Sc_ARI)
mean(Pois_P_SBM_ARI)
mean(Pois_B_SBM_ARI)
# performance similar to Poisson SBM when Poisson data. our model doesn't use multiple 
# restarts, etc
# poisson sbm works well with hawkes data.

### Compare how fast online learns ####

sbm_count = BM_poisson('SBM',A_matrix, explore_min = K, explore_max = K, verbosity = 0)
sbm_count$estimate()
sbm_count_est = apply(sbm_count$memberships[[3]]$Z,1,which.max)
return(adjustedRandIndex(true_Z,sbm_count_est))

# compare results from 
early = results_online$early_tau
final = results_online$tau
early_est = apply(early,1,which.max)
final_est = apply(final, 1, which.max)
adjustedRandIndex(early_est,final_est)
adjustedRandIndex(Z,final_est)
adjustedRandIndex(Z,early_est)
