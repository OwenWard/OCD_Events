#### Simulation demonstrating online recovery of
# communitites ####

library(Rcpp)
library(RcppArmadillo)
library(mclust)
sourceCpp("onlineblock.cpp")

set.seed(200)

n_sims = 50
# store online nmi for each of the simulations 

Time <- 500
dT <- 5
inter_T <- 5
inters <- Time/(dT*inter_T)
sim_nmi = matrix(nrow=n_sims,ncol=inters)


for(iter in 1:n_sims) {
  print(iter)
  Time = 500
  #dT = 10
  K <- 3
  Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
  B_Pois <- matrix(0,K,K,byrow = TRUE)
  m <- 100
  Pi <- matrix(c(0.4,0.3,0.3),1,3)
  Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))
  
  A <- list()
  for(i in 1:m){
    # could sample these here with SBM structure...
    edge_list <- c(1:m)
    edge <- sample(edge_list[-i], 40) - 1
    edge <- sort(edge)
    A[[i]] <- edge
  }
  
  system.time(alltimes <- sampleBlockHak(Time, A, Z, Mu, B_Pois, lam = 1))
  
  
  
  Pi = rep(1/K,K)
  B = matrix(runif(K*K),K,K)
  diag(B) = rnorm(3,mean = 1, sd = 0.1)
  tau = matrix(runif(m*K),nrow=m,ncol=K)
  tau = tau/rowSums(tau)
  S = matrix(0,nrow = m,ncol = K)
  #num_wind = Time/dT
  #inter_num = num_wind/50
  #dT <- 10
  #inter_T <- 1
  
  #results_online <- estimate_Poisson(full_data = alltimes,tau,B,Pi,S,A,m,K,dT=0.1,Time)
  results_online <- estimate_Poisson(full_data = alltimes,
                                     A,m,K,Time,dT,B,
                                     tau,Pi,S,inter_T)
  print(results_online$Pi)
  time_points = dim(results_online$early_tau)[3]-1
  early_nmi = rep(0,time_points)
  for(i in 1:time_points){
    est_Z = apply(results_online$early_tau[,,i],1,which.max)
    early_nmi[i]= aricode::NMI(Z,est_Z)
  }
  
  #online_nmi = cbind(windows = seq(50,500,50),early_nmi)
  online_nmi = cbind(iters = seq(from = dT*inter_T,
                                 by = dT*inter_T,
                                 length.out = time_points),early_nmi)
  sim_nmi[iter,] = early_nmi
}



#as_tibble(online_nmi) %>% ggplot(aes(iters,early_nmi)) + geom_line()

apply(sim_nmi,2,mean) 

# then want to plot this somehow
library(tidyverse)

# nmi_online = tibble(windows = seq(50,500,50),nmi = t(sim_nmi))
# 
# nmi_online$windows
# nmi_online




new_sim_nmi = t(sim_nmi)
colnames(new_sim_nmi) = paste("Sim",c(1:n_sims),sep="")



new_sim_nmi

nmi_online = as_tibble(new_sim_nmi )
nmi_online %>% head()
nmi_online$Windows = seq(dT*inter_T,500,dT*inter_T)

nmi_online %>% head()

# then pivot longer?

nmi_online %>% pivot_longer(
  cols = starts_with("Sim")) %>%
  ggplot(aes(Windows,value)) +
  geom_point() + geom_smooth() + ylab("NMI")


nmi_online %>% pivot_longer(
  cols = starts_with("Sim")) %>%
  group_by(Windows) %>%
  summarise(Ave = mean(value),sdev = sd(value)) %>%
  mutate(se = sdev/sqrt(50)) %>%
  ggplot(aes(Windows,Ave)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=Ave-se, ymax=Ave+se), width=.1,color="blue") +
  ylab("Average NMI") + xlab("Time") +
  theme(axis.title = element_text(size = 14))
saveRDS(nmi_online,"Simulations/nmi_online.RDS")


# analysis of results from Sim 1 also ###
sim1_t = readRDS("Simulations/output_sim1_t.RDS")
sim1_t

sim1_t$nmi

new_nmi_t = t(sim1_t$nmi)
n_sims = 50
colnames(new_nmi_t) = paste("Sim",c(1:n_sims),sep="")
t_nmi = as_tibble(new_nmi_t)
t_nmi$t = sim1_t$Time_setting
t_nmi %>% pivot_longer(
  cols = starts_with("Sim")) %>%
  group_by(t) %>%
  summarise(Ave = mean(value),sdev = sd(value)) %>%
  mutate(se = sdev/sqrt(50))  %>%
  ggplot(aes(t,Ave)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin = Ave-se,ymax = Ave + se),width=.1, color="blue") +
  ylab("Average NMI") + xlab("Observation Time")



####
sim1_m = readRDS("Simulations/output_sim1_m.RDS")
# only did 25 simulations in the end
nmi_m = sim1_m$nmi[1:25,] 
nmi_m

new_nmi_m = t(nmi_m)
n_sims = 25
colnames(new_nmi_m) = paste("Sim",c(1:n_sims),sep="")
m_nmi = as_tibble(new_nmi_m)
m_nmi$m = sim1_m$m

m_nmi %>% pivot_longer(
  cols = starts_with("Sim")) %>%
  group_by(m) %>%
  summarise(Ave = mean(value),sdev = sd(value)) %>%
  mutate(se = sdev/sqrt(25)) %>%
  ggplot(aes(m,Ave)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=Ave-se, ymax=Ave+se), width=.1,color="blue") +
  ylab("Average NMI") + xlab("Number Nodes")
  
