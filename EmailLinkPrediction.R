library(Rcpp)
library(RcppArmadillo)
sourceCpp("C:/Users/owenw/Documents/Online_Point_Process/onlineblock.cpp")



library(tidyverse)


# --- ####
emails <- readRDS("emailscleaned.RDS")

emails_train = emails %>% filter(Time < 471)
summary(emails_train$Time)
summary(emails$Time)

emails_test = emails %>% filter(Time > 471)

A = list()

for(i in 1:nrow(emails)){
  j = emails$Send[i]
  k = emails$Rec[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

for(i in 1:nrow(emails)){
  j = emails$Send[i]
  k = emails$Rec[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}


A_test = lapply(A,unique)
#A_test = A_test[-1]

max_time = max(emails$Time)
max_time/(3600*24) # so in seconds, this transformation puts it in
# days

### need to recode ####
#summary(as.numeric(factor(emails$Send,levels = users)))
m = length(unique(c(emails$Send,emails$Rec)))

# then fit the poisson model 
# problem with the indexing of the users

Time = 804
dT = 0.5
K = 3
Pi = c(0.3,0.3,0.4)
B = matrix(c(2,0.5,0.5,0.5,1.1,.65,0.75,0.85,1.15),nrow = K,ncol = K,byrow = T)
B = matrix(runif(K*K),K,K)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)





# Homogeneous Poisson ####
train_time = 400
test_time = 500
results_train <- estimate_Poisson(full_data = as.matrix(emails_train),tau,
                                  B,Pi,S,A_test,m,K,dT=1.5,T=train_time)

# then predict links
est_Z = apply(results_train$tau,1,which.max)


# simulate data and get counts #
B_Pois = matrix(0,nrow = K,ncol = K)
pred_times = sampleBlockHak_pre(T = test_time,startT = train_time,A_test,est_Z-1,
                                Mu = results_train$B,B_Pois,lam=1)

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)


emails_train = emails %>% filter(Time < train_time)
summary(emails_train$Time)
summary(emails$Time)

emails_test = emails %>% filter(Time > train_time) %>%
  filter(Time < test_time)

summary(emails_test$Time)

pred_set = pred_times %>% group_by(Send,Rec) %>% tally() %>%
  filter(Send != Rec)
test_set = emails_test %>% group_by(Send,Rec) %>% tally()

### simulating events
test_set %>% left_join(pred_set,by = c("Send","Rec")) %>%
  mutate(n.y = replace_na(n.y,0)) %>%
  rename(Pred_Count_Sim = n.y,
         num_test = n.x) %>%
  mutate(diff_sim = Pred_Count_Sim-num_test,
         diff_zer = 0-num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_sim^2)),
            RMSE_0 = sqrt(mean(diff_zer^2))) 

# approx 15-30
# this is comparable to the models is Miscourdiou et al

# compute mean number of events also 
groups = tibble(Send = c(1:m)-1,Z = est_Z)

pred_test_set =  test_set %>% left_join(groups, by = "Send")

groups = tibble(Rec = c(1:m)-1,Z = est_Z)
pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")



mu = rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] = results_train$B[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set$Mean = mu


pred_test_set %>% 
  mutate(mean_events = Mean*(test_time-train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))
# this gives comparable results, but interesting that the null model does
# so well


## Homogeneous Hawkes Model ####
K = 4
dT = 3
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results_hawkes_sim <- online_estimator_eff(as.matrix(emails_train), 
                                           A_test, m, K, T = train_time, dT, lam = 1, B, Mu, tau)
# 
itermax <- 802 / dT
stop_eps <- 0.001
results.batch <- batch_estimator(as.matrix(emails), A_test, m, K, T=802, dT, 
                                 lam = 1.0, B, Mu, tau, itermax, stop_eps)
#
results_hawkes_sim = results.batch

est_Z = apply(results_hawkes_sim$tau,1,which.max)

pred_times = sampleBlockHak_pre(T = test_time,startT = train_time,A_test,est_Z-1,
                                Mu = results_hawkes_sim$Mu,
                                B = results_hawkes_sim$B,
                                lam = results_hawkes_sim$lam)

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)

pred_set = pred_times %>% group_by(Send,Rec) %>% tally() %>% 
  filter(Send!=Rec)



#full_prediction using same as neurips 18 paper
test_set %>% left_join(pred_set,by = c("Send","Rec")) %>%
  mutate(n.y = replace_na(n.y,0)) %>%
  rename(Pred_Count_Sim = n.y,
         num_test = n.x) %>%
  mutate(diff_sim = Pred_Count_Sim-num_test,
         diff_zer = 0-num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_sim^2)),
            RMSE_0 = sqrt(mean(diff_zer^2))) 
# these are marginally different but not much

##### using mean number of events ####
# new prediction function


est = Predict_Counts_Hawkes(startT = 472,finalT = 804,A_test,est_Z-1,
                            results_hawkes_sim$Mu,results_hawkes_sim$B,m,
                            results_hawkes_sim$lam)

est = est[-1,]
colnames(est)= c("Send","Rec","Time")
est = as_tibble(est)

test_set %>% left_join(est,by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))
