.libPaths("/rigel/stats/users/ogw2103/rpackages")
setwd("/rigel/stats/users/ogw2103/code/Online_Point_Process")
library(Rcpp)
library(RcppArmadillo)
library(mclust)
sourceCpp("onlineblock.cpp")


num_Sims = 5
Times = c(25,50,100,200)
m_setting = c(50,100,200,500)
K = 3


### Varying m with fixed T = 100
results = list(adj = matrix(,nrow=num_Sims,ncol = 4),
               nmi = matrix(,nrow = 50,ncol = 4),
               m = m_setting,
               Time_setting = Times,
               B = array(dim = c(K,K,num_Sims,4)))



for(setting_m in c(1:4)){
  #for(setting_t in c(1:4)){
    Time = 100#Times[setting_t]
    cat("Time Length: ",Time, "\n")
    m = m_setting[setting_m]
    cat("Num Nodes: ",m, "\n")
    cat("Num Sims: ", num_Sims, "\n")
    for(iter in c(1:num_Sims)){
      #print(iter)
      #print(Time)
      dT = 0.5
      K <- 3
      Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
      B_Pois <- matrix(0,K,K,byrow = TRUE)
      Pi <- matrix(c(0.4,0.3,0.3),1,3)
      Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))
      
      A <- list()
      for(i in 1:m){
        # could sample these here with SBM structure...
        num_edge = m*0.4
        edge <- sample(m, num_edge) - 1
        edge <- sort(edge)
        A[[i]] <- edge
      }
      
      alltimes <- sampleBlockHak(Time, A, Z, Mu, B_Pois, lam = 1)
      Pi = rep(1/K,K)
      B = matrix(runif(K*K),K,K)
      #diag(B) = rnorm(3,mean = 1, sd = 0.1)
      tau = matrix(runif(m*K),nrow=m,ncol=K)
      tau = tau/rowSums(tau)
      S = matrix(0,nrow = m,ncol = K)
      
      
      results_online <- estimate_Poisson(full_data = alltimes,
                                         tau,B,Pi,S,A,m,K,dT=0.1,Time)
      
      #result = norm(Mu-results_online$B)
      est_Z = apply(results_online$tau,1,which.max)
      ari = adjustedRandIndex(Z,est_Z)
      nmi = aricode::NMI(Z,est_Z)
      results$B[,,iter,setting_m] = results_online$B
      results$adj[iter,setting_m] = ari
      results$nmi[iter,setting_m] = nmi
    } 
  #}
}


saveRDS(results,"output_sim1_m.RDS")


### Varying m with fixed T = 100
results2 = list(adj = matrix(,nrow=num_Sims,ncol = 4),
               nmi = matrix(,nrow = 50,ncol = 4),
               m = m_setting,
               Time_setting = Times,
               B = array(dim = c(K,K,num_Sims,4)))




for(setting_t in c(1:4)){
  Time = 100#Times[setting_t]
  cat("Time Length: ",Time, "\n")
  m = m_setting[setting_m]
  cat("Num Nodes: ",m, "\n")
  cat("Num Sims: ", num_Sims, "\n")
  for(iter in c(1:num_Sims)){
    #print(iter)
    #print(Time)
    dT = 0.5
    K <- 3
    Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
    B_Pois <- matrix(0,K,K,byrow = TRUE)
    Pi <- matrix(c(0.4,0.3,0.3),1,3)
    Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))
    
    A <- list()
    for(i in 1:m){
      # could sample these here with SBM structure...
      num_edge = m*0.4
      edge <- sample(m, num_edge) - 1
      edge <- sort(edge)
      A[[i]] <- edge
    }
    
    alltimes <- sampleBlockHak(Time, A, Z, Mu, B_Pois, lam = 1)
    Pi = rep(1/K,K)
    B = matrix(runif(K*K),K,K)
    #diag(B) = rnorm(3,mean = 1, sd = 0.1)
    tau = matrix(runif(m*K),nrow=m,ncol=K)
    tau = tau/rowSums(tau)
    S = matrix(0,nrow = m,ncol = K)
    
    
    results_online <- estimate_Poisson(full_data = alltimes,
                                       tau,B,Pi,S,A,m,K,dT=0.5,Time)
    
    #result = norm(Mu-results_online$B)
    est_Z = apply(results_online$tau,1,which.max)
    ari = adjustedRandIndex(Z,est_Z)
    nmi = aricode::NMI(Z,est_Z)
    results2$B[,,iter,setting_m] = results_online$B
    results2$adj[iter,setting_m] = ari
    results2$nmi[iter,setting_m] = nmi
  } 
}


saveRDS(results2,"output_sim1_t.RDS")