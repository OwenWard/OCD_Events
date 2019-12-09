library(Rcpp)
library(RcppArmadillo)
sourceCpp("C:/Users/owenw/Documents/Online_Point_Process/onlineblock.cpp")


setwd("C:/Users/owenw/Downloads/")

library(dplyr)
emails = read.csv(gzfile("email-Eu-core-temporal.txt.gz"))

emails %>% head()

colnames(emails) = c("Data")

#emails %>% separate(sep = " ")

#separate(emails,Data," ")

emails = emails %>% separate(Data,c("Send","Rec","Time"),sep = " ")


emails = emails %>% mutate(Time = as.numeric(Time))
emails = emails %>% arrange(Time)




# then fit the models on both and compare link predictions
users = unique(c(emails$Send,emails$Rec))
m = length(unique(c(emails$Send,emails$Rec)))
m

# need to construct the corresponding A for this data


emails = emails %>% mutate(Send = as.numeric(Send),
                           Rec = as.numeric(Rec))

emails = emails %>% 
  mutate(Send = as.numeric(factor(Send,levels = users))-1)
emails = emails %>% 
  mutate(Rec = as.numeric(factor(Rec,levels = users))-1)


emails = emails %>% mutate(Time = Time/(3600*24))

# ---
emails <- readRDS("emailscleaned.RDS")

emails_train = emails %>% filter(Time < 471)
summary(emails_train$Time)
summary(emails$Time)

emails_test = emails %>% filter(Time > 471)

A = list()

for(i in 1:m){
  A[[i+1]] = i
}

for(i in 1:nrow(emails)){
  j = emails$Send[i]
  k = emails$Rec[i]
  #print(j)
  #A[[j+1]] = j
  #A[[k+1]] = k
  A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}
A_test = lapply(A,unique)
#A_test = A_test[-1]

max_time = max(emails$Time)
max_time/(3600*24) # so in seconds, this transformation puts it in
# days

### need to recode ####
summary(as.numeric(factor(emails$Send,levels = users)))


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



results <- estimate_Poisson(full_data = as.matrix(emails),tau,
                            B,Pi,S,A_test,m,K,dT,T=Time)
results$B
plot(x=c(1:length(results$AveELBO)),y=results$AveELBO,type = 'l')

# try do link prediction
results_train <- estimate_Poisson(full_data = as.matrix(emails_train),tau,
                            B,Pi,S,A_test,m,K,dT=0.5,T=471)

# then predict links
est_Z = apply(results_train$tau,1,which.max)


# not the real predicted times...
B_Pois = matrix(0,nrow = K,ncol = K)
pred_times = sampleBlockHak_pre(T = Time,startT = 471,A,est_Z-1,
                                Mu = results_train$B,B_Pois,lam=1)

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)

pred_set = pred_times %>% group_by(Send,Rec) %>% tally()
test_set = emails_test %>% group_by(Send,Rec) %>% tally()


comparison = test_set %>% full_join(pred_set,by = c("Send","Rec") )

comparison = comparison %>% mutate(n.x = replace_na(n.x,0),
                                   n.y = replace_na(n.y,0),
                                   diff = n.x-n.y)

comparison$n.x[is.na(comparison$n.x)] = 0
comparison$n.y[is.na(comparison$n.y)] = 0

sqrt(mean(comparison$diff^2))
# approx 500
# this is comparable to the models is Miscourdiou et al


## fit the Hom Hawkes Model also
K = 2
dT = 5
Pi = rep(1/K,K)
B = matrix(runif(K*K),K,K)
Mu = matrix(runif(K*K),K,K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau = matrix(runif(m*K),nrow=m,ncol=K)
tau = tau/rowSums(tau)
S = matrix(0,nrow = m,ncol = K)


system.time(results_hawkes_sim <- online_estimator_eff(as.matrix(emails_train),
                                                       A_test, m, K, T = 471, dT, 
                                                       lam = 1, B, Mu, tau))


dT <- 2.25
K <- 2
b <- 0.5
W1 <- matrix(0,m,K)
W2 <- matrix(0,m,K)
system.time(results.ccrm <- ccrm_estimator(alltimes,A,m,K,T,dT,lam = 1.0,W1,W2,b))
# seems to blow up for certain dT?

est_Z = apply(results_hawkes_sim$tau,1,which.max)

pred_times = sampleBlockHak_pre(T = Time,startT = 472,A,est_Z-1,
                                Mu = results_hawkes_sim$Mu,
                                B = results_hawkes_sim$B,
                                lam = results_hawkes_sim$lam)

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)

pred_set = pred_times %>% group_by(Send,Rec) %>% tally()
test_set = emails_test %>% group_by(Send,Rec) %>% tally()


comparison = test_set %>% full_join(pred_set,by = c("Send","Rec") )

comparison = comparison %>% mutate(n.x = replace_na(n.x,0),
                                   n.y = replace_na(n.y,0),
                                   diff = n.x-n.y)

comparison$n.x[is.na(comparison$n.x)] = 0
comparison$n.y[is.na(comparison$n.y)] = 0
comparison$diff <- comparison$n.x - comparison$n.y


sqrt(mean(comparison$diff^2))
# approx 61.2791 for k = 3

## Non Homogeneous Poisson ####
K = 2
H = 1
window = 1

results_npois_sim <- nonhomoPois_estimator(as.matrix(emails_train),A_test,m,K,H,
                                            window,T=471,dT=2.5, gravity = 0.0,MuA,tau)
est_Z = apply(results_nhawkes_sim$tau,1,which.max)
# then predictions for this

pred_times = sampleBlockHak_nonhomo_pre(T = Time, startT = 472,A,est_Z-1,
                                        MuA = results_nhawkes_sim$MuA,
                                        B = results_nhawkes_sim$B,
                                        window,
                                        lam = 1)





## Non Homogeneous Hawkes ####
K = 2
H <- 2
window = 1/H

MuA = array(0,dim=c(K,K,H))

results_nhawkes_sim <- nonhomoHak_estimator_eff(as.matrix(emails_train),A_test,m,K,H,
                                       window,T=471,dT=2.5,lam = 0.1, gravity = 0.0, B,MuA,tau)
est_Z = apply(results_nhawkes_sim$tau,1,which.max)

pred_times = sampleBlockHak_nonhomo_pre(T = Time, startT = 472,A,est_Z-1,
                                        MuA = results_nhawkes_sim$MuA,
                                        B = results_nhawkes_sim$B,
                                        window,
                                        lam = results_nhawkes_sim$lam)

colnames(pred_times) = c("Send","Rec","Time")
pred_times = as_tibble(pred_times)

pred_set = pred_times %>% group_by(Send,Rec) %>% tally()
test_set = emails_test %>% group_by(Send,Rec) %>% tally()


comparison = test_set %>% full_join(pred_set,by = c("Send","Rec") )

comparison = comparison %>% mutate(n.x = replace_na(n.x,0),
                                   n.y = replace_na(n.y,0),
                                   diff = n.x-n.y)

comparison$n.x[is.na(comparison$n.x)] = 0
comparison$n.y[is.na(comparison$n.y)] = 0
comparison$diff <- comparison$n.x - comparison$n.y

sqrt(mean(comparison$diff^2))

