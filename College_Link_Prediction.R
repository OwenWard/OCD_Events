# this is a tidy file which consists of link prediction
# metric for any of the datasets on the snap website, given that
# they follow the standard format


library(Rcpp)
library(RcppArmadillo)
sourceCpp("C:/Users/owenw/Documents/Online_Point_Process/onlineblock.cpp")

### Read in raw data and preprocess to our format ####
# this section uses dplyr
library(tidyverse)
#college = read.csv(gzfile("C:/Users/owenw/Downloads/CollegeMsg.txt.gz"))
college = read.csv(gzfile("C:/Users/owenw/Downloads/sx-mathoverflow.txt.gz"))
# although we call the dataset college here in can be any dataset


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
saveRDS(college,"math_overflow_messages.RDS")


#### analyse the data ####
college = readRDS("math_overflow_messages.RDS")
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
Time = 2351
K = 2
dT = 10
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
#train_time = 76.38
train_time = 1934
#test_time = 193
test_time = 2351

#### Hom Poisson ####
dT = 50
K = 2
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)

results_pois_train <- estimate_Poisson(full_data = as.matrix(college_train),tau,
                                 B,Pi,S,A_test,m,K,dT,T=train_time)
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
K = 3
dT = 25
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


results_hawkes_sim <- online_estimator_eff(as.matrix(college_train), 
                                           A_test, m, K, T = train_time, dT, lam = 1, B, Mu, tau)

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
 