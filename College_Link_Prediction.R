# this is a tidy file which consists of link prediction
# metric for any of the datasets on the snap website, given that
# they follow the standard format


library(Rcpp)
library(RcppArmadillo)
sourceCpp("C:/Users/owenw/Documents/Online_Point_Process/onlineblock.cpp")

### Read in raw data and preprocess to our format ####
# this section uses dplyr
library(tidyverse)
college = read.csv(gzfile("C:/Users/owenw/Downloads/CollegeMsg.txt.gz"))
#college = read.csv(gzfile("C:/Users/owenw/Downloads/sx-mathoverflow.txt.gz"))
# although we call the dataset college here in can be any dataset
#college = read.csv(gzfile("C:/Users/owenw/Downloads/email-Eu-core-temporal.txt.gz"))


dim(college)

colnames(college) = c("Data")

#emails %>% separate(sep = " ")

#separate(emails,Data," ")

college = college %>% separate(Data,c("Send","Rec","Time"),sep = " ")

college %>% head()

college = college %>% mutate(Time = as.numeric(Time))
college = college %>% arrange(Time)

summary(college$Time/(3600*24))
summary(college$Time)
range(college$Time)
(max(college$Time) - min(college$Time))/(3600*24)
# so this is coded in days, just with an offset term at the start

college = college %>% mutate(Time = Time - min(Time))
college = college %>% mutate(Time = Time/(3600*24))
# now range from 0 to 193 days


# then transform user so it is coded up in useable format
users = unique(c(college$Send,college$Rec))
m = length(unique(c(college$Send,college$Rec)))
m # same as on website

college = college %>% 
  mutate(Send = as.numeric(factor(Send,levels = users))-1) %>%
  mutate(Rec = as.numeric(factor(Rec,levels = users))-1)

summary(college$Send)
# this looks good


hist(college$Time)
# this is all over the place in terms of email distribution also
saveRDS(college,"college_messages.RDS")


#### analyse the data ####
#college = readRDS("math_overflow_messages.RDS")
A = list()

m = length(unique(c(college$Send,college$Rec)))

for(i in 1:nrow(college)){
  j = college$Send[i]
  k = college$Rec[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

for(i in 1:nrow(college)){
  j = college$Send[i]
  k = college$Rec[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]],college$Rec[i])
}


A_test = lapply(A,unique)


# check it fits on all data 
Time = 804
K = 2
dT = 2
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

results_pois <- estimate_Poisson(full_data = as.matrix(college),tau,
                            B,Pi,S,A_test,m,K,dT,T=Time)

results_hawkes_sim <- online_estimator_eff(as.matrix(college), 
                                           A_test, m, K, T = Time, dT, lam = 1, B, Mu, tau)

results_hawkes_sim$Pi

# split into training and test

college_train = college %>% top_frac(-.85)
# this gives a training time cutoff of 1934 for the big dataset
# total time is 2350
college_test = college %>% top_frac(.15)
summary(college_train$Time)
summary(college_test$Time)
summary(college$Time)

# then do link prediction on the test set
#train_time = 470.3
train_time = 76.38
#train_time = 1934
#test_time = 804
test_time = 193
#test_time = 2351

#### Hom Poisson ####
dT = 0.5 # such that approx 400 windows for whole time period
#dT = .25
K = 2 # 2 for college, 4 for email
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

# online estimate
results_pois_train <- estimate_Poisson(full_data = as.matrix(college_train),tau,
                                 B,Pi,S,A_test,m,K,dT,T=train_time)

# batch estimate for hom-Poisson needed here


# simulate events
est_Z = apply(results_pois_train$tau,1,which.max)
B_Pois = matrix(0,nrow = K,ncol = K)
pred_times = sampleBlockHak_pre(T = test_time,startT = train_time,
                                A_test,est_Z-1,
                                Mu = results_pois_train$B,B_Pois,lam=1)

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)
pred_set = pred_times %>% group_by(Send,Rec) %>% tally() %>%
  filter(Send != Rec)
test_set = college_test %>% group_by(Send,Rec) %>% tally()

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


# mean events
groups = tibble(Send = c(1:m)-1,Z = est_Z)

pred_test_set =  test_set %>% left_join(groups, by = "Send")

groups = tibble(Rec = c(1:m)-1,Z = est_Z)
pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")



mu = rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] = results_pois_train$B[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
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


### Hom Hawkes ####
K = 2
dT = 0.5
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

# online estimator
results_hawkes_sim <- online_estimator_eff(as.matrix(college_train), 
                                           A_test, m, K, T = train_time, dT, lam = 1, B, Mu, tau)

# batch estimator..
results_hawkes_sim <- batch_estimator(as.matrix(college_train), 
                                           A_test, m, K, T = train_time, dT, lam = 1, B, Mu, tau,
                                      itermax =  1000,stop_eps = 0.01)

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

# mean number of events
# this is very slow for the large network. will try make it faster
est = Predict_Counts_Hawkes(startT = train_time,finalT = test_time,A_test,est_Z-1,
                            results_hawkes_sim$Mu,results_hawkes_sim$B,m,
                            results_hawkes_sim$lam)

est = est[-1,]
colnames(est)= c("Send","Rec","Time")
est = as_tibble(est)

test_set %>% left_join(est,by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))
 

## InHomogeneous Poisson ####

# taking the mean of the time component?
window = 1/7
K <- 3
H <- 7
dT = 10
MuA_start = array(0,c(K,K,H))
tau_start = matrix(0,m,K)
non_hom_pois_est <- nonhomoPois_estimator(as.matrix(college_train),A_test,m,K,H,window,
                                          T = train_time,dT,gravity = 0.001,MuA_start,tau_start)

# taking the average of these basis functions for link prediction
dim(non_hom_pois_est$MuA)
baseline = apply(non_hom_pois_est$MuA,c(1,2),mean)

est_Z = apply(non_hom_pois_est$tau,1,which.max)
est_Z

groups = tibble(Send = c(1:m)-1,Z = est_Z)

pred_test_set =  test_set %>% left_join(groups, by = "Send")

groups = tibble(Rec = c(1:m)-1,Z = est_Z)
pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")



mu = rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] = baseline[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set$Mean = mu

pred_test_set %>% 
  mutate(mean_events = Mean*(test_time-train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))


### InHomogeneous Hawkes ####


B_start = matrix(0,K,K)
non_homo_Hawkes_est = nonhomoHak_estimator_eff(as.matrix(college_train),A_test,m,K,H,window, lam = 1,
                                               T = train_time,dT,gravity = 0.001,MuA_start = MuA_start,
                                               tau_start = tau_start,B_start = B_start)
non_homo_Hawkes_est$MuA
non_homo_Hawkes_est$Pi
non_homo_Hawkes_est$B


baseline = apply(non_homo_Hawkes_est$MuA,c(1,2),mean)

est_Z = apply(non_homo_Hawkes_est$tau,1,which.max)
est_Z

est = Predict_Counts_Hawkes(startT = train_time,finalT = test_time,A_test,est_Z-1,
                            baseline,non_homo_Hawkes_est$B,m,
                            non_homo_Hawkes_est$lam)

est = est[-1,]
colnames(est)= c("Send","Rec","Time")
est = as_tibble(est)

test_set %>% left_join(est,by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))

### Asymmetric CCRM Predictions ####

W1 <- matrix(0,m,K)
W1[,1] <- c(rep(1,m))
W2 <- matrix(0,m,K)
W2[,1] <- c(rep(1,m))
b <- 0.5

results.ccrm_train <- ccrm_estimator(as.matrix(college_train),A_test,m,K,T=train_time,
                               dT,lam = 1.0,W1,W2,b)

#dim(results.ccrm_train$W1)
# now have individual baseline for each node pair...
mu = rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] = results.ccrm_train$W1[pred_test_set$Send[i]+1,]%*%results.ccrm_train$W2[pred_test_set$Rec[i]+1,]
}

min(mu)
k = results.ccrm_train$lam*(1-results.ccrm_train$b)
exp_term = exp(-k*train_time) - exp(-k*test_time)
b = results.ccrm_train$b
lam = results.ccrm_train$lam
mean_events = (test_time-train_time)-b/(lam*(1-b)^2)*exp_term
mean_events_ccrm = mean_events*mu/(1-b)


pred_test_set %>% 
  add_column(mean_events_ccrm = mean_events_ccrm) %>%
  #mutate(mean_events = mean_events_ccrm) %>%
  mutate(diff_mean = n -mean_events_ccrm) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

### CCRM debug on training data...

results.ccrm_train$W1[16+1,]%*%results.ccrm_train$W2[283+1,]
results.ccrm_train$W1[180,]


as_tibble(results.ccrm_train$W1) %>% mutate(max_row = rowSums) %>% arrange(desc(max_row))

rowSums(resu)


