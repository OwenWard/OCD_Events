library(Rcpp)
library(RcppArmadillo)
library(mclust)
sourceCpp("onlineblock.cpp")

set.seed(200)


n_sims <- 50
sim_nmi = matrix(nrow=n_sims,ncol=20)


for(iter in 1:n_sims){
  Time = 200
  dT = 1.25
  K <- 3
  Mu_H <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.2,0.75),K,K,byrow = TRUE)
  B_H <- matrix(c(0.5,0.1,0.3,0.4,0.4,0.4,0.2,0.6,0.2),K,K,byrow = TRUE)
  m <- 100
  Pi <- matrix(c(0.4,0.3,0.3),1,3)
  Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))
  
  A <- list()
  for(i in 1:m){
    edge <- sample(m, 25) - 1
    edge <- sort(edge)
    A[[i]] <- edge
  }
  
  system.time(alltimes_hawkes <- sampleBlockHak(Time, A, Z, Mu_H,
                                                B_H, lam = 1))
  
  
  Pi = rep(1/K,K)
  Mu = matrix(runif(K*K),K,K)
  B = matrix(runif(K*K),K,K)
  #diag(B) = rnorm(3,mean = 1, sd = 0.1)
  tau = matrix(runif(m*K),nrow=m,ncol=K)
  tau = tau/rowSums(tau)
  S = matrix(0,nrow = m,ncol = K)
  
  
  results_hawkes_sim <- online_estimator_eff_revised(alltimes_hawkes, 
                                                     A, m, K, Time, dT=1,
                                                     lam = 1, B, Mu, tau,
                                                     inter_T = 10, 
                                                     is_elbo = F)
  time_points = dim(results_hawkes_sim$early_tau)[3]-1
  early_nmi = rep(0,time_points)
  for(i in 1:time_points){
    est_Z = apply(results_hawkes_sim$early_tau[,,i],1,which.max)
    early_nmi[i]= aricode::NMI(Z,est_Z)
  }
  
  sim_nmi[iter,] = early_nmi
  
  
}


# then analyse this online NMI
apply(sim_nmi,2,mean) 

# then want to plot this somehow
library(tidyverse)

# nmi_online = tibble(windows = seq(50,500,50),nmi = t(sim_nmi))
# 
# nmi_online$windows
# nmi_online




new_sim_nmi = t(sim_nmi)
colnames(new_sim_nmi) = paste("Sim",c(1:n_sims),sep="")

dT = 1
inter_T = 10
Time = 200

new_sim_nmi

nmi_online = as_tibble(new_sim_nmi )
nmi_online %>% head()
nmi_online$Windows = seq(inter_T,200,inter_T)

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


