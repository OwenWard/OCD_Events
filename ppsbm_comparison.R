generated_Q3_n20$z


results = mainVEM(data = generated_Q3_n20$data, n = 20, Qmin = 3, Qmax = 3, method = 'kernel',
                  directed = F,
                  d_part = 1, n_perturb = 1)
est_Z = apply(results[[1]]$tau, 2, which.max)
true_Z = apply(generated_Q3_n20$z, 2, which.max)
adjustedRandIndex(true_Z,est_Z)

dim(results[[1]]$logintensities.ql.ij)
#### simulate new data ####
prop.groups <- c(0.5,0.5)

# 3 different intensity functions :
intens <- list(NULL)
intens[[1]] <- list(intens= function(x) 100*x*exp(-8*x),max=5)
# (q,l) = (1,1)
intens[[2]] <- list(intens= function(x) exp(3*x)*(sin(6*pi*x-pi/2)+1)/2,max=13)
# (q,l) = (1,2)
intens[[3]] <- list(intens= function(x) 8.1*(exp(-6*abs(x-1/2))-.049),max=8)
# (q,l) = (2,2)
intens[[4]] <- list(intens= function(x) 1.1*(exp(-2*abs(x-1/2))-.1),max=3)

obs <- generateDynppsbm(intens,Time=1,n=20,prop.groups,directed=T)


### fit this with our model also. reformat for this ####
# need to use this and then match across
# difficulty is this data being undirected
n = 20
index = listNodePairs(n,directed = T)
Id = convertNodePair(index[,1],index[,2],n,directed = T)
Id 
# so these correspond exactly in some sense

n_events = length(obs$data$time.seq)

#generated_Q3_n20$data$type.seq

alltimes_ppsbm = matrix(0,nrow = n_events,ncol = 3)
for(i in 1:n_events){
  ind = obs$data$type.seq[i]
  alltimes_ppsbm[i,1] = index[ind,1]-1
  alltimes_ppsbm[i,2] = index[ind,2]-1
  alltimes_ppsbm[i,3] = obs$data$time.seq[i]
}

alltimes_ppsbm # can then use this
# need to construct the corresponding A also
A_ppsbm = matrix(0,nrow = n,ncol = n)

#A_ppsbm = list()
for(i in 1:n_events){
  j = alltimes_ppsbm[i,1] # sender id
  k = alltimes_ppsbm[i,2] # rec id
  # print(j)
  # print("sends to")
  # print(k)
  A_ppsbm[j,k] = A_ppsbm[j,k] + 1
  # A_ppsbm[[j]] = c(A_ppsbm[[j]],j-1)
  # A_ppsbm[[k]] = c(A_ppsbm[[k]],k-1)
  # A_ppsbm[[j]] = c(A_ppsbm[[j]],k-1)
}

# this A is not correct
A_pp = list()
for(i in 1:n){
  temp = which(A_ppsbm[i,] > 0)-1
  A_pp[[i]] = c(i-1,temp)
}


#####
K = 2
m = 20
T = 1
dT = 0.01
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
diag(B) = rnorm(2,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

Mu = matrix(runif(K*K),K,K)


results_online <- estimate_Poisson(full_data = alltimes_ppsbm,tau,B,Pi,S,A_pp,m,K,dT,T)
results_hawkes <- online_estimator(alltimes_ppsbm, A_pp, m, K, T, dT, lam = 1, B, Mu, tau)

true_Z = apply(obs$z,2,which.max)

est_Z_online_Poisson = apply(results_online$tau,1,which.max)
adjustedRandIndex(est_Z_online_Poisson,true_Z)


est_Z_online_H = apply(results_hawkes$tau,1,which.max)
adjustedRandIndex(est_Z_online_H,true_Z)


results = mainVEM(data = obs$data, n = 20, Qmin = 2, Qmax = 2, method = 'kernel',
                  directed = T,
                  d_part = 1, n_perturb = 1)
est_Z = apply(results[[1]]$tau, 2, which.max)
adjustedRandIndex(true_Z,est_Z)
