# because these models don't scale, use on a subset of the data?
library(mclust)
library(dplyr)
library(ppsbm)

library(Rcpp)
library(RcppArmadillo)
sourceCpp("onlineblock.cpp")


#### Simulate data from our Model, fit with both methods ####
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


results_online <- estimate_Poisson(full_data = alltimes,tau,B,Pi,S,A,m,K,dT,T)
results_online$B
Mu
est_Z = apply(results_online$tau,1,which.max)#-1
est_Z
Z


adjustedRandIndex(Z,est_Z)

# then transform and fit with PPSBM
alltimes %>% head()

n_events = dim(alltimes)[1]
ids = rep(0,n_events)
times = rep(0,n_events)

for(i in 1:n_events){
  j = alltimes[i,1]+1
  k = alltimes[i,2]+1
  ids[i] = convertNodePair(j,k,m,TRUE)
  times[i] = alltimes[3]
}

sim_data = list(time.seq = times, type.seq = ids, Time = T)

results = mainVEM(data = sim_data, n = m, Qmin = 3, Qmax = 3, method = 'kernel',
                  d_part = 1, n_perturb = 1)

results[[2]]$tau
est_Z_ppsbm = apply(results[[1]]$tau, 2, which.max)
est_Z_ppsbm
adjustedRandIndex(Z,est_Z_ppsbm)
# so with


##
##
##### simulate from ppsbm ####
# simulate data from this model
generated_Q3

convertGroupPair(c(1:3),c(1:3),Q=3)

generated_Q3_n20$z
generated_Q3_n20$data

n_events = length(generated_Q3_n20$data$type.seq)

for(i in 1:n_events){
  convertNodePair(1,2,n=20)
  convertGroupPair(3,3,3)
}


# Generate data from an undirected graph with n=10 individuals and Q=2 clusters

# equal cluster proportions
prop.groups <- c(0.5,0.5)

# 3 different intensity functions :
intens <- list(NULL)
intens[[1]] <- list(intens= function(x) 100*x*exp(-8*x),max=5)
# (q,l) = (1,1)
intens[[2]] <- list(intens= function(x) exp(3*x)*(sin(6*pi*x-pi/2)+1)/2,max=13)
# (q,l) = (1,2)
intens[[3]] <- list(intens= function(x) 8.1*(exp(-6*abs(x-1/2))-.049),max=8)
# (q,l) = (2,2)

# generate data :
obs <- generateDynppsbm(intens,Time=10,n=30,prop.groups,directed=FALSE)

obs$data$time.seq

results = mainVEM(data = obs$data, n = 100, Qmin = 2, Qmax = 5, method = 'kernel')


generated_Q3_n20
results = mainVEM(data = generated_Q3_n20$data,n=20,Qmin=2,Qmax = 3,
                  method = 'kernel',n_perturb = 1,d_part = 1)
results