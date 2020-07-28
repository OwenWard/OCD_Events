#### Simulation to show the difficulties of binning data ####

### code for cluster
#.libPaths("/moto/stats/users/ogw2103/rpackages/")
#setwd("/moto/stats/users/ogw2103/Code/OCD_Events/Sims/")
n_sims <- 100
###


library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

setwd("C:/Users/owenw/Documents/Github_Local/OCD_Events/")
sourceCpp("cpp_files/onlineblock.cpp")
#sourceCpp("../cpp_files/onlineblock.cpp")

source("Sims/pensky_fcns.R")

bin_fun <- function(events,m,max_Time,window_size){
  winds <- seq(from=window_size,to=max_Time,by= window_size)
  adj <- array(0, dim = c(m,m,length(winds)))
  for(i in 1:length(winds)){
    # get edges in each window and then updating the corresponding
    # entry in the adjacency matrix
    upper <- winds[i]
    lower <- upper-window_size
    wind_events <- events %>% 
      filter( Time < upper & Time > lower) %>%
      group_by(start,end) %>% tally()
    for(j in 1:nrow(wind_events)){
      row <- wind_events$start[j]+1
      col <- wind_events$end[j] + 1
      adj[row,col,i] <- wind_events$n[j] 
    }
  }
  return(adj)
}

theme_set(theme_classic())

set.seed(200)

# ##### Inhomogeneous Hawkes ####

###n_sims <- 100
m_values <- c(50,100,250,500)
no_methods <- 3 #ours, sum all, pz_estimator


l <- 2
m <- 2 # need m =l for final time point
m0 <- 1
l0 <- 2

# results <- matrix(NA, nrow = n_sims,ncol = no_methods+2 )
#
# all_results <- list(m_50 = results, m_100 = results, m_200 = results,
#                     m_500 = results)
all_results <- list(rep(tibble(),length(m_values)))
#all_results <- tibble()

for(i in seq_along(m_values)) {
  m_nodes <- m_values[i]
  #results <- matrix(NA, nrow = n_sims,ncol = no_methods )
  # then each sim
  results <- matrix(NA, nrow = n_sims,ncol = no_methods+2 )
  colnames(results) <- c("SC_Sum","OCD","SC_PZ","nodes","sim")
  results <- as_tibble(results)
  for(j in 1:n_sims) {
    cat("SIM ",j,"\n",sep = "")
    K <- 2
    H <- 2
    MuA <- array(0,c(K,K,H))
    MuA[,,1] <- matrix(c(0.8,0.2,0.6,0.4),2,2)
    MuA[,,2] <- matrix(c(0.4,0.7,0.2,0.7),2,2)

    Mu <- matrix(c(0.6,0.2,0.3,0.5),K,K,byrow = TRUE)
    B <- matrix(c(0.1,0.02,0.01,0.05),K,K,byrow = TRUE)
    Pi <- matrix(c(0.4,0.6),1,2)
    Z <- c(rep(0,m_nodes*Pi[1]),rep(1,m_nodes*Pi[2]))
    Time <- 100

    A <- list()
    for(k in 1:m_nodes){
      # could sample these here with SBM structure...
      node_list <- c(1:m_nodes)
      edge <- sample(node_list[-k], 10) - 1
      edge <- sort(edge)
      A[[k]] <- edge
    }

    system.time(alltimes <- sampleBlockHak_nonhomo(Time,
                                                   A, Z, MuA, B,
                                                   window = 1, lam = 1))
    ### fit the models
    B_start <- matrix(runif(K*K),K,K)
    MuA_start <- array(runif(K*K*H),c(K,K,H))

    tau_start = matrix(runif(m_nodes*K),nrow=m_nodes,ncol=K)
    tau_start = tau_start/rowSums(tau_start)

    nonHawkes_online <- nonhomoHak_estimator_eff_revised(alltimes,
                                                         A,m_nodes,K,H,
                                                         window = 0.75,
                                                         T = Time,dT=10,
                                                         lam = 1.75,
                                                         gravity = 0.001,
                                                         B_start,
                                                         MuA_start,
                                                         tau_start)

    est_Z <- apply(nonHawkes_online$tau,1,which.max)


    test_events <- tibble(start = alltimes[,1],end = alltimes[,2],
                          Time = alltimes[,3])

    out <- bin_fun(test_events,m_nodes,max_Time = Time,window_size = 1)

    ### store the results
    # using sum of all events

    #results$sum[j] <- aricode::NMI(fcd::spectral.clustering(out[,,dim(out)[3]]),Z)
    # using only events in final time window
    sum_adj <- apply(out, c(1,2), sum)#/dim(out)[3]
    results$SC_Sum[j] <- aricode::NMI(fcd::spectral.clustering(sum_adj,K=2),Z)
    # pz
    pz_est <- pz_estimator_3(out,time = Time,
                             l0 = l0,
                             m0 = m0,
                             m = m,
                             r= Time-1)
    results$SC_PZ[j] <- aricode::NMI(fcd::spectral.clustering(pz_est,K=2),Z)
    # our method
    results$OCD[j] <- aricode::NMI(Z,est_Z)
    results$nodes[j] <- m_nodes
    results$sim[j] <- j

  }
  all_results[[i]] <- results

}


saveRDS(all_results, file = "Sims/output_sims/agg_comp_inhom_hawk_wind_85.RDS")

##### Plot Inhom Hawkes ####
plot_results <- all_results %>%
  bind_rows() %>%
  pivot_longer(cols=SC_Sum:SC_PZ, names_to = "Method") %>%
  group_by(nodes,Method) %>%
  summarise(NMI = mean(value), se = sd(value)/sqrt(n_sims)) %>%
  rowwise() %>%
  mutate(lower = max(NMI - se,0), upper = min(NMI + se,1)) %>%
  ggplot(aes(nodes,NMI,colour = Method)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper))

plot_results

#### Block Hawkes ####


####

####n_sims <- 20
m_values <- c(50,100,250,500)
no_methods <- 3 #ours, sum all, pz_estimator


l <- 2
m <- 2 # need m =l for final time point
m0 <- 1
l0 <- 2


all_results <- list(rep(tibble(),length(m_values)))
#all_results <- tibble()

for(i in seq_along(m_values)) {
  m_nodes <- m_values[i]
  #results <- matrix(NA, nrow = n_sims,ncol = no_methods )
  # then each sim
  results <- matrix(NA, nrow = n_sims,ncol = no_methods+2 )
  colnames(results) <- c("SC_Sum","OCD","SC_PZ","nodes","sim")
  results <- as_tibble(results)
  for(j in 1:n_sims) {
    cat("SIM ",j," M ",m_nodes, "\n",sep = "")
    K <- 2
    # H <- 2
    # MuA <- array(0,c(K,K,H))
    # MuA[,,1] <- matrix(c(0.8,0.2,0.6,0.4),2,2)
    # MuA[,,2] <- matrix(c(0.4,0.7,0.2,0.7),2,2)

    Mu <- matrix(c(0.6,0.2,0.3,0.5),K,K,byrow = TRUE)
    B <- matrix(c(0.1,0.02,0.01,0.05),K,K,byrow = TRUE)
    Pi <- matrix(c(0.4,0.6),1,2)
    Z <- c(rep(0,m_nodes*Pi[1]),rep(1,m_nodes*Pi[2]))
    Time <- 200

    A <- list()
    for(k in 1:m_nodes){
      # could sample these here with SBM structure...
      node_list <- c(1:m_nodes)
      edge <- sample(node_list[-k], 10) - 1
      edge <- sort(edge)
      A[[k]] <- edge
    }

    system.time(alltimes <- sampleBlockHak(Time,A, Z, Mu, B, lam = 1))
    ### fit the models
    B_start <- matrix(runif(K*K),K,K)
    # MuA_start <- array(runif(K*K*H),c(K,K,H))
    Mu_start <- matrix(runif(K*K),K,K)

    tau_start = matrix(runif(m_nodes*K),nrow=m_nodes,ncol=K)
    tau_start = tau_start/rowSums(tau_start)

    Hawkes_online <- online_estimator_eff_revised(alltimes,
                                                         A,m_nodes,K,
                                                         T = Time,dT=5,
                                                         lam = 1.25,
                                                         B_start,
                                                         Mu_start,
                                                         tau_start,inter_T = 1)

    est_Z <- apply(Hawkes_online$tau,1,which.max)


    test_events <- tibble(start = alltimes[,1],end = alltimes[,2],
                          Time = alltimes[,3])

    out <- bin_fun(test_events,m_nodes,max_Time = Time,window_size = 1)

    ### store the results
    # using sum of all events

    #results$sum[j] <- aricode::NMI(fcd::spectral.clustering(out[,,dim(out)[3]]),Z)
    # using only events in final time window
    sum_adj <- apply(out, c(1,2), sum)#/dim(out)[3]
    results$SC_Sum[j] <- aricode::NMI(fcd::spectral.clustering(sum_adj,K=2),Z)
    # pz
    pz_est <- pz_estimator_3(out,time = Time,
                             l0 = l0,
                             m0 = m0,
                             m = m,
                             r= Time-1)
    results$SC_PZ[j] <- aricode::NMI(fcd::spectral.clustering(pz_est,K=2),Z)
    # our method
    results$OCD[j] <- aricode::NMI(Z,est_Z)
    results$nodes[j] <- m_nodes
    results$sim[j] <- j

  }
  all_results[[i]] <- results

}

all_results %>%
  bind_rows() %>%
  pivot_longer(cols=SC_Sum:SC_PZ, names_to = "Method") %>%
  group_by(nodes,Method) %>%
  summarise(NMI = mean(value), se = sd(value)/sqrt(n_sims)) %>%
  rowwise() %>%
  mutate(lower = max(NMI - se,0), upper = min(NMI + se,1)) %>%
  ggplot(aes(nodes,NMI,colour = Method)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper))

saveRDS(all_results, file = "Sims/output_sims/agg_comp_hawk.RDS")

#### Block Poisson ####

####n_sims <- 5
m_values <- c(50,100,250,500)
no_methods <- 3 #ours, sum all, pz_estimator


l <- 2
m <- 2 # need m =l for final time point
m0 <- 1
l0 <- 2


all_results <- list(rep(tibble(),length(m_values)))
#all_results <- tibble()

for(i in seq_along(m_values)) {
  m_nodes <- m_values[i]
  #results <- matrix(NA, nrow = n_sims,ncol = no_methods )
  # then each sim
  results <- matrix(NA, nrow = n_sims,ncol = no_methods+2 )
  colnames(results) <- c("SC_Sum","OCD","SC_PZ","nodes","sim")
  results <- as_tibble(results)
  for(j in 1:n_sims) {
    print(j)
    K <- 2
    # H <- 2
    # MuA <- array(0,c(K,K,H))
    # MuA[,,1] <- matrix(c(0.8,0.2,0.6,0.4),2,2)
    # MuA[,,2] <- matrix(c(0.4,0.7,0.2,0.7),2,2)
    
    Mu <- matrix(c(0.6,0.2,0.3,0.5),K,K,byrow = TRUE)
    B <- matrix(0,K,K,byrow = TRUE)
    Pi <- matrix(c(0.4,0.6),1,2)
    Z <- c(rep(0,m_nodes*Pi[1]),rep(1,m_nodes*Pi[2]))
    Time <- 100
    
    A <- list()
    for(k in 1:m_nodes){
      # could sample these here with SBM structure...
      node_list <- c(1:m_nodes)
      edge <- sample(node_list[-k], 10) - 1
      edge <- sort(edge)
      A[[k]] <- edge
    }
    
    system.time(alltimes <- sampleBlockHak(Time,A, Z, Mu, B, lam = 1))
    ### fit the models
    Pi = rep(1/K,K)
    B = matrix(runif(K*K),K,K)
    #diag(B) = rnorm(K,mean = 1, sd = 0.1)
    tau = matrix(runif(m_nodes*K),nrow=m_nodes,ncol=K)
    tau = tau/rowSums(tau)
    S = matrix(0,nrow = m_nodes,ncol = K)
    
    poisson_online <- estimate_Poisson(full_data = alltimes,
                     A,m_nodes,K,Time,dT = 5,
                     step_size = 0.5,
                     B,
                     tau,Pi,S,inter_T = 1,is_elbo = FALSE)
    
    est_Z <- apply(poisson_online$tau,1,which.max)
    
    
    test_events <- tibble(start = alltimes[,1],end = alltimes[,2],
                          Time = alltimes[,3])
    
    out <- bin_fun(test_events,m_nodes,max_Time = Time,window_size = 1)
    
    ### store the results
    # using sum of all events
    
    #results$sum[j] <- aricode::NMI(fcd::spectral.clustering(out[,,dim(out)[3]]),Z) 
    # using only events in final time window
    sum_adj <- apply(out, c(1,2), sum)#/dim(out)[3]
    results$SC_Sum[j] <- aricode::NMI(fcd::spectral.clustering(sum_adj,K=2),Z)
    # pz
    pz_est <- pz_estimator_3(out,time = Time,
                             l0 = l0,
                             m0 = m0,
                             m = m,
                             r= Time-1)
    results$SC_PZ[j] <- aricode::NMI(fcd::spectral.clustering(pz_est,K=2),Z)
    # our method
    results$OCD[j] <- aricode::NMI(Z,est_Z)
    results$nodes[j] <- m_nodes
    results$sim[j] <- j
    
  }
  all_results[[i]] <- results
  
}


saveRDS(all_results, file = "Sims/output_sims/agg_comp_pois.RDS")
beepr::beep()

# all_results %>%
#   bind_rows() %>%
#   pivot_longer(cols=SC_Sum:SC_PZ, names_to = "Method") %>%
#   group_by(nodes,Method) %>%
#   summarise(NMI = mean(value), se = sd(value)/sqrt(n_sims)) %>%
#   rowwise() %>%
#   mutate(lower = max(NMI - se,0), upper = min(NMI + se,1)) %>%
#   ggplot(aes(nodes,NMI,colour = Method)) + 
#   geom_line() +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) + 
#   scale_x_continuous(breaks = m_values)

