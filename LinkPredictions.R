### Link Prediction

source("comparison_fcns.R")


#
#### Simulate data from our Model, fit with all methods ####
T = 50
dT = 0.1
K <- 3
Mu <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)

B_Pois <- matrix(0,K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  edge <- sample(m, 40) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak(T, A, Z, Mu_H, B_H, lam = 1))
#new_times = sampleBlockHak_pre(T = 20,startT = 10,A,Z,Mu,B_Pois,lam=1)

colnames(alltimes) = c("Send","Rec","Time")
alltimes = as_tibble(alltimes)


train_set = alltimes %>% filter(Time < 40)
test_set = alltimes %>% filter(Time > 40 )

# no issue down to here

Pi = c(0.3,0.3,0.4)
B = matrix(c(1,0.5,0.5,0.5,1,.65,0.75,0.85,1.),nrow = K,ncol = K,byrow = T)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

# more random start
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results_train <- estimate_Poisson(full_data = as.matrix(train_set),tau,B,Pi,S,A,m,K,dT,T=40)
results_train$B

results_hawkes_sim <- online_estimator(as.matrix(train_set), A, m, K, Time = 40, dT=0.2, lam = 1, B, Mu, tau)

Mu
est_Z = apply(results_train$tau,1,which.max)
adjustedRandIndex(Z,est_Z)

## or down to here....
### 
### has to be here then

# not the real predicted times...
pred_times = sampleBlockHak_pre(T = 50,startT = 40,A,est_Z-1,
                               Mu = results_train$B,B_Pois,lam=1)


# and ok down to here...
colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)

saveRDS(pred_times,"pred_events.RDS")
saveRDS(test_set,"test_events.RDS")
#### down to here is fine ....

test_counts = test_set %>% group_by(Send,Rec) %>% tally()
pred_counts = pred_times %>% gr


# reloading and not having run above 
test_set = readRDS("true_events.RDS")
pred_set = readRDS("pred_events.RDS")

library(dplyr)
library(tidyverse)

pred_set = pred_set %>% group_by(Send,Rec) %>% tally()
test_set = test_set %>% group_by(Send,Rec) %>% tally()


comparison = test_set %>% full_join(pred_set,by = c("Send","Rec") )

#comparison = 
  
comparison = comparison %>% mutate(n.x = replace_na(n.x,0),
                                   n.y = replace_na(n.y,0),
                                   diff = n.x-n.y)


sqrt(mean(comparison$diff^2))



#### Homogeneous Hawkes Predictions ####

Time = 500
train_time = 400
final_time = 500

dT = 0.1
K <- 3

Mu_H <- matrix(c(0.6,0.2,0.3,0.1,1.0,0.4,0.5,0.4,0.8),K,K,byrow = TRUE)
B_H <- matrix(c(0.8,0.2,0.1,0.15,0.75,0.25,0.1,0.25,0.9),K,K,byrow = TRUE)

m <- 500
Pi <- matrix(c(0.4,0.3,0.3),1,3)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]),rep(2,m*Pi[3]))

A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  edge <- sample(m, 200) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


alltimes_hawkes <- sampleBlockHak(Time, A, Z, Mu_H, B_H, lam = 1)

colnames(alltimes_hawkes) = c("Send","Rec","Time")
alltimes_hawkes <- as_tibble(alltimes_hawkes)


test_set = alltimes_hawkes %>% filter(Time > train_time)
train_set = alltimes_hawkes %>% filter(Time < train_time)

# then fit on training data

Pi = rep(1/K,K)
Mu = matrix(runif(K*K),K,K)
B = matrix(runif(K*K),K,K)
diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results_hawkes_test <- online_estimator(as.matrix(train_set), A, m, K, T=train_time, dT=0.2, lam = 1, B, Mu, tau)

# estimate from the full data
itermax <- Time / dT
stop_eps <- 0.005
results_hawkes_batch <- batch_estimator(as.matrix(train_set), A, m, K, T= train_time, dT, 
                          lam = 1.0, B, Mu, tau, itermax, stop_eps)


est_Z = apply(results_hawkes_test$tau,1,which.max)-1
est_Z_full = apply(results_hawkes_batch$tau,1,which.max)-1

pred_times = sampleBlockHak_pre(T=final_time,startT = train_time,A,est_Z,
                                Mu = results_hawkes_test$Mu,
                                B = results_hawkes_test$B,
                                lam = results_hawkes_test$lam)

pred_times_batch = sampleBlockHak_pre(T=final_time,startT = train_time,A,est_Z_full,
                                Mu = results_hawkes_batch$Mu,
                                B = results_hawkes_batch$B,
                                lam = results_hawkes_batch$lam)


saveRDS(test_set,"true_events.RDS")
saveRDS(pred_times,"pred_events.RDS")
saveRDS(pred_times_batch,"pred_events_full.RDS")

library(tidyverse)

test_set = readRDS("true_events.RDS")
pred_times = readRDS("pred_events.RDS")
pred_times_batch = readRDS("pred_events_full.RDS")

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)

colnames(pred_times_batch) = c("Send","Rec","Time")
pred_times_batch = as_tibble(pred_times_batch)

colnames(test_set) = c("Send","Rec","Time")
test_set = as_tibble(test_set)

test_counts = test_set %>% group_by(Send,Rec) %>% tally()
pred_counts = pred_times %>% group_by(Send,Rec) %>% tally()
pred_counts_batch = pred_times_batch %>% group_by(Send,Rec) %>% tally()


comparison = test_counts %>% full_join(pred_counts,by = c("Send","Rec") )
comparison2 = test_counts %>% full_join(pred_counts_batch,
                                        by=c("Send","Rec"))
#comparison = 

comparison = comparison %>% mutate(n.x = replace_na(n.x,0),
                                   n.y = replace_na(n.y,0),
                                   diff = n.y-n.x)

comparison2 = comparison2 %>% mutate(n.x = replace_na(n.x,0),
                                   n.y = replace_na(n.y,0),
                                   diff = n.y-n.x)

sqrt(mean(comparison$diff^2))
sqrt(mean(comparison2$diff^2))


### Non Homogeneous Hawkes

##

#### Hawkes Non-Homogeneous ####
Time <- 50
train_time = 40
K <- 2
H <- 2
MuA <- array(0,c(K,K,H))
MuA[,,1] <- matrix(c(0.8,0.2,0.6,0.4),2,2)
MuA[,,2] <- matrix(c(0.4,0.7,0.2,0.7),2,2)
B <- matrix(c(0.8,0.2,0.4,0.6),K,K,byrow = TRUE)
m <- 100
Pi <- matrix(c(0.6,0.4),1,K)
Z <- c(rep(0,m*Pi[1]),rep(1,m*Pi[2]))
window <- 0.25

A <- list()
for(i in 1:m){
  edge <- sample(m, 20) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}



system.time(alltimes <- sampleBlockHak_nonhomo(T=Time, A, Z, MuA, B, window, lam = 1))

dT <- 2.5
tau <- matrix(0,m,K)
for (k in 1:K){
  tau[which(Z == (k-1)),k] <- 1
}

## training and test set
colnames(alltimes) = c("Send","Rec","Time")
alltimes = as_tibble(alltimes)

train_set = alltimes %>% filter(Time < train_time)
test_set = alltimes %>% filter(Time > train_time)




system.time(results.online <- nonhomoHak_estimator(as.matrix(train_set),A,m,K,H,
                                                   window,T=Time,dT,lam = 0.1, gravity = 0.0, B,MuA,tau))

itermax <- T / dT
stop_eps <- 0.001
system.time(results.batch <- batch_nonhomoHak_estimator(as.matrix(train_set),A,m,K,H,
                                                        window,T=Time,dT,lam = 0.1, gravity = 0.0,
                                                        B,MuA,tau,itermax, stop_eps))

est_Z = apply(results.online$tau,1,which.max)-1
est_Z_full = apply(results.batch$tau,1,which.max)-1
