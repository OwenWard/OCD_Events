#### Oral Link Prediction Example ####
# would like to add a plot showing how the model 
# improves in terms of link prediction over time...




# break the data up into 3 hour windows, fit, predict for the
# window, then add that data and fit
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(lubridate)
sourceCpp("C:/Users/owenw/Desktop/Online_Point_Process/onlineblock.cpp")
#set.seed(200)



nyc_data <- read_csv("C:/Users/owenw/Google Drive/Tian/Oral/data/201907-citibike-tripdata.csv.zip")

# can use starttime or stoptime here
citi_data <- nyc_data %>% dplyr::select(`start station id`,`end station id`,stoptime,`start station latitude`,
                                        `start station longitude`,
                                        `end station latitude`,
                                        `end station longitude`) %>% arrange(stoptime)

citi_data %>% ggplot(aes(stoptime)) + geom_histogram()

rm(nyc_data)

summary(citi_data$stoptime)
# so data over approximately one month
# want to convert this into days
start_time <- ymd_hms("2019-07-01 00:00:00", tz = "UTC")
citi_data$minutes <- citi_data$stoptime-start_time

head(citi_data$minutes)

#citi_data_full <- citi_data # so don't have to read in every time

# this will actually be hours now
citi_data <- citi_data %>% mutate(days = as.numeric(minutes/60)) %>%
  filter(days < 24) # for the full month 24*31



# sometimes this minutes is in minutes and sometimes this
# is in seconds.... 
# not sure how you can specify it. check with
dim(citi_data)


citi_data %>% group_by(`start station id`,`end station id`) %>% 
  tally() %>% arrange(desc(n)) %>% 
  filter(`start station id`!=`end station id`) 

#citi_data %>% ggplot(aes(days)) + geom_histogram()
#citi_data %>% ggplot(aes(stoptime)) + geom_histogram()

# so a pretty small proportion of these trips are between the same
# stations

citi_data <- citi_data  %>%
  filter(`start station id`!= `end station id`) %>%
  filter(`end station id`!= 3791) # jersey city

bike_data <- citi_data %>% drop_na() %>% 
  dplyr::select(`start station id`,`end station id`,days)
# then convert the stations to be numeric factors which start at 0


bike_data$`start station id` <- as.numeric(as.factor(bike_data$`start station id`))-1
bike_data$`end station id` <- as.numeric(as.factor(bike_data$`end station id`))-1

bike_data %>% group_by(`start station id`,`end station id`) %>% tally()

m <- length(unique(c(bike_data$`start station id`,
                     bike_data$`end station id`)))

citi_data$start_stat <- bike_data$`start station id` +1
citi_data$end_stat <- bike_data$`end station id`+1
# can use this to join them after





A = list()

for(i in 1:nrow(bike_data)){
  j = bike_data$`start station id`[i]
  k = bike_data$`end station id`[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
  #A[[j+1]] = c(A[[j+1]],emails$Rec[i])
}

#init_A <- A

for(i in 1:nrow(bike_data)){
  j = bike_data$`start station id`[i]
  k = bike_data$`end station id`[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]],bike_data$`end station id`[i])
}




A_test = lapply(A,unique)

rm(A)


### do this for homogeneous Poisson initially ####
n_sims <- 1
time_breaks <- seq(from=3,to=24,by=3)

inter_T <- 1

output <- tibble()

for(sims in 1:n_sims){
  print(sims)
  for(j in 1:(length(time_breaks)-1)){
    print(j)
    end_train = time_breaks[j]
    train_data <- bike_data %>% filter(days < end_train)
    end_test = time_breaks[j+1]
    test_data <- bike_data %>% filter(days > end_train & days < end_test)
    K <- 4
    Time = end_train
    dT <- 0.1
    Pi = rep(1/K,K)
    B = matrix(runif(K*K),K,K)
    tau = matrix(runif(m*K),nrow=m,ncol=K)
    tau = tau/rowSums(tau)
    S = matrix(0,nrow = m,ncol = K)
    online_pois <- estimate_Poisson(full_data = as.matrix(train_data),
                                    #full_data = as.matrix(sample_data),
                                    A_test,m,K,
                                    T=Time,
                                    dT,B,
                                    tau,Pi,S,inter_T,is_elbo = TRUE)
    test_set <- test_data %>% 
      rename(Send=`start station id`,Rec = `end station id`)%>%
      group_by(Send,Rec) %>% tally()
    est_Z <- apply(online_pois$tau,1,which.max)
    
    groups = tibble(Send = c(1:m)-1,Z = est_Z)
    
    pred_test_set =  test_set %>% left_join(groups, by = "Send")
    
    groups = tibble(Rec = c(1:m)-1,Z = est_Z)
    pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")
    
    mu = rep(0,dim(pred_test_set)[1])
    for(i in 1:length(mu)){
      mu[i] = online_pois$B[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
    }
    
    pred_test_set$Mean = mu
    output_current <- pred_test_set %>% 
      mutate(mean_events = Mean*(end_test-end_train)) %>%
      mutate(diff_mean = n -mean_events) %>%
      mutate(diff_zero = n - 0) %>%
      ungroup() %>%
      summarise(RMSE = sqrt(mean(diff_mean^2)),
                RMSE_0 = sqrt(mean(diff_zero^2)))
    print(output_current)
    output_current$end_train <- end_train
    output_current$sim <- sims
      output <- output %>% bind_rows(output_current)
  }
  
}

# then plot this....
output

output %>% 
  pivot_longer(-c(end_train,sim),names_to = "Model",values_to = "RMSE") %>%
  ggplot(aes(end_train,RMSE,color=Model)) + geom_smooth(se=FALSE) 


output %>% 
  mutate(diff = RMSE_0-RMSE) %>% ggplot(aes(end_train,diff)) + 
  geom_boxplot(aes(group=end_train))



### repeat for inhomogeneous Poisson ####
### could even just do this for homogeneous Poisson initially
n_sims <- 2
time_breaks <- seq(from=3,to=24,by=3)

output_nh_pois <- tibble()

for(sims in 1:n_sims){
  print(sims)
  for(j in 1:(length(time_breaks)-1)){
    print(j)
    end_train = time_breaks[j]
    train_data <- bike_data %>% filter(days < end_train)
    end_test = time_breaks[j+1]
    test_data <- bike_data %>% filter(days > end_train & days < end_test)
    K <- 4
    Time = end_train
    dT <- 0.1
    window = 1/4
    H <- 4
    dT = 0.25 # 2 for email, 0.5 for college, 6 for math
    MuA_start = array(0,c(K,K,H))
    tau_start = matrix(0,m,K)
    non_hom_pois_est <- nonhomoPois_estimator(as.matrix(train_data),A_test,m,K,H,window,
                                              T = Time,dT,gravity = 0.001,
                                              MuA_start,tau_start)
    
    test_set <- test_data %>% 
      rename(Send=`start station id`,Rec = `end station id`)%>%
      group_by(Send,Rec) %>% tally()
    est_Z <- apply(non_hom_pois_est$tau,1,which.max)
    baseline = apply(non_hom_pois_est$MuA,c(1,2),mean)
    
    groups = tibble(Send = c(1:m)-1,Z = est_Z)
    
    pred_test_set =  test_set %>% left_join(groups, by = "Send")
    
    groups = tibble(Rec = c(1:m)-1,Z = est_Z)
    pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")
    
    mu = rep(0,dim(pred_test_set)[1])
    for(i in 1:length(mu)){
      mu[i] = baseline[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
    }
    
    pred_test_set$Mean = mu
    output_current <- pred_test_set %>% 
      mutate(mean_events = Mean*(end_test-end_train)) %>%
      mutate(diff_mean = n -mean_events) %>%
      mutate(diff_zero = n - 0) %>%
      ungroup() %>%
      summarise(RMSE = sqrt(mean(diff_mean^2)),
                RMSE_0 = sqrt(mean(diff_zero^2))) 
    print(output_current)
    output_current$end_train <- end_train
    output_current$sim <- sims
    output_nh_pois <- output_nh_pois %>% bind_rows(output_current)
  }
  
}

# then plot this....
output_nh_pois

output_nh_pois %>% 
  pivot_longer(-c(end_train,sim),names_to = "Model",values_to = "RMSE") %>%
  ggplot(aes(end_train,RMSE,color=Model)) + geom_smooth(se=FALSE) 


output_nh_pois %>% 
  mutate(diff = RMSE_0-RMSE) %>% ggplot(aes(end_train,diff)) + 
  geom_boxplot(aes(group=end_train))

# do I want to plot homog vs nhomog against each other on the one plot,
# not show the RMSE_0 at all?



#### Predict out degree from each station 
### Homogeneous Poisson ####
n_sims <- 1
wind <- 3
time_breaks <- seq(from=wind,to=24,by=wind)


output <- tibble()
output_edge <- tibble()



for(sims in 1:n_sims){
  print(sims)
  for(j in 1:(length(time_breaks)-1)){
    print(j)
    end_train = time_breaks[j]
    start_train = end_train-wind
    train_data <- bike_data %>% 
      filter(days < end_train & days > start_train)
    train_data_batch <- bike_data %>% 
      filter(days < end_train)
    end_test = time_breaks[j+1]
    test_data <- bike_data %>% filter(days > end_train & days < end_test)
    K <- 5
    Time = end_train
    dT <- 0.1
    if(j == 1){
      Pi = rep(1/K,K)
      B = matrix(runif(K*K),K,K)
      tau = matrix(runif(m*K),nrow=m,ncol=K)
      tau = tau/rowSums(tau)
      S = matrix(0,nrow = m,ncol = K)
    }
    if(j!= 1){
      B = online_pois$B
      tau = online_pois$tau
      Pi = online_pois$Pi
      S = online_pois$S
    }

    online_pois <- estimate_Poisson(full_data = as.matrix(train_data),
                                    #full_data = as.matrix(sample_data),
                                    A_test,m,K,
                                    T=Time,
                                    dT,B,
                                    tau,Pi,S,inter_T,is_elbo = TRUE)
    test_set <- test_data %>% 
      rename(Send=`start station id`,Rec = `end station id`)%>%
      group_by(Send,Rec) %>% tally()
    # for out degree
    test_set_send <- test_set %>% group_by(Send) %>% tally()
    
    est_Z <- apply(online_pois$tau,1,which.max)
    
    groups = tibble(Send = c(1:m)-1,Z = est_Z)
    
    pred_test_set =  test_set %>% left_join(groups, by = "Send")
    
    groups = tibble(Rec = c(1:m)-1,Z = est_Z)
    pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")
    
    mu_online = rep(0,dim(pred_test_set)[1])
    for(i in 1:length(mu_online)){
      mu_online[i] = online_pois$B[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
    }
    pred_test_set$Mean_online = mu_online
    # output_current <- pred_test_set %>% 
    #   mutate(mean_events = Mean*(end_test-end_train)) %>%
    #   mutate(diff_mean = n -mean_events) %>%
    #   mutate(diff_zero = n - 0) %>%
    #   ungroup() %>%
    #   summarise(RMSE = sqrt(mean(diff_mean^2)),
    #             RMSE_0 = sqrt(mean(diff_zero^2)))
    
    
    
    
    #-------------------------------#
    # repeat for batch estimator #
    B_batch = matrix(runif(K*K),K,K)
    tau_batch = matrix(runif(m*K),nrow=m,ncol=K)
    tau_batch = tau/rowSums(tau)
    S_batch = matrix(0,nrow = m,ncol = K)
    batch_pois <-batch_estimator_hom_Poisson(as.matrix(train_data_batch),
                                             A_test,m,K,
                                             T=Time,
                                             B_batch,tau_batch,               
                                             itermax = 50,
                                             stop_eps = 0.01)
    est_Z_batch <- apply(batch_pois$tau,1,which.max)
    
    groups = tibble(Send = c(1:m)-1,Z = est_Z_batch)
    pred_test_set =  pred_test_set %>% left_join(groups, by = "Send")
    
    groups = tibble(Rec = c(1:m)-1,Z = est_Z_batch)
    pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")
    mu_batch = rep(0,dim(pred_test_set)[1])
    for(i in 1:length(mu_batch)){
      mu_batch[i] = batch_pois$B[pred_test_set$Z.x.x[i],pred_test_set$Z.y.y[i]]
    }
    pred_test_set$Mean_batch = mu_batch
    
    
    
    # -------- combine the results from both ---- #
    pred_test_set_events <- pred_test_set %>% 
      mutate(mean_events_online = Mean_online*(end_test-end_train),
             mean_events_batch = Mean_batch*(end_test-end_train))
    
    
    out_pred <- pred_test_set_events %>% group_by(Send) %>%
      mutate(true_out = sum(n),est_online = sum(mean_events_online),
             est_batch = sum(mean_events_batch)) %>%
      dplyr::select(Send,true_out,est_online,est_batch) %>%
      distinct(,.keep_all = T) %>%
      ungroup() %>%
      summarise(mean_resid_online = 
                  mean( (true_out-est_online)/sqrt(est_online) ),
                mean_resid_batch = 
                  mean( (true_out-est_batch)/sqrt(est_batch) )
                )
      # summarise(RMSE_online = sqrt(mean((true_out-est_online)^2)),
      #           RMSE_batch = sqrt(mean((true_out-est_batch)^2)))
    
    ### want a better metric for this if we can
    
    # print(output_current)
    # output_current$end_train <- end_train
    # output_current$sim <- sims
    # output <- output %>% bind_rows(output_current)
    out_pred$end_train <- end_train
    out_pred$sim <- sims
    output_edge <- output_edge %>% bind_rows(out_pred)
    
  }
  
}

beepr::beep()
# now would need to pivot_longer to compare online and batch here
#output_edge %>% ggplot(aes(end_train,RMSE)) + geom_smooth()
#output_edge$end_train

output_edge %>% 
  pivot_longer(-c(end_train,sim),names_to = "Model",values_to = "RMSE") %>%
  ggplot(aes(end_train,RMSE,color=Model)) + geom_smooth() +
  xlab("Hours history") + 
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22)) + 
  scale_color_hue(labels=c("batch","online"))+
  #theme_bw() + 
  theme(axis.title = element_text(size=14))

################################
# ============================
# ============================
# same procedure for inhomogeneous Poisson


n_sims <- 20
wind <- 1.5
time_breaks <- seq(from=wind,to=24,by=wind)
H <- 4


output_nh <- tibble()
output_edge_nh <- tibble()



for(sims in 1:n_sims){
  print(sims)
  for(j in 1:(length(time_breaks)-1)){
    print(j)
    end_train = time_breaks[j]
    start_train = end_train-wind
    train_data <- bike_data %>% 
      filter(days < end_train & days > start_train)
    train_data_batch <- bike_data %>% 
      filter(days < end_train)
    end_test = time_breaks[j+1]
    test_data <- bike_data %>% filter(days > end_train & days < end_test)
    K <- 4
    Time = end_train
    dT <- 0.1
    if(j == 1){
      MuA_start = array(runif(K*K*H),c(K,K,H))
      tau_start = matrix(runif(m*K),m,K) #
      tau_start = t(apply(tau_start,1,function(x)(x/sum(x))))
      Pi_start = rep(1/K,K)
      S_start = matrix(0,nrow = m,ncol = K)
    }
    if(j!= 1){
      MuA_start = non_hom_pois_est$MuA
      tau_start = non_hom_pois_est$tau
      Pi_start = non_hom_pois_est$Pi
      S_start = non_hom_pois_est$S
    }
    
    window = 1/4
    dT = 0.25 
    
    # MuA_start = array(0,c(K,K,H))
    # tau_start = matrix(0,m,K)
    
    
    
    non_hom_pois_est <- nonhomoPois_estimator_update(as.matrix(train_data),A_test,m,K,H,window,
                                              T = Time,dT,gravity = 0.001,
                                              MuA_start,tau_start,
                                              S_start)
    test_set <- test_data %>% 
      rename(Send=`start station id`,Rec = `end station id`)%>%
      group_by(Send,Rec) %>% tally()
    # for out degree
    test_set_send <- test_set %>% group_by(Send) %>% tally()
    
    est_Z <- apply(non_hom_pois_est$tau,1,which.max)
    
    groups = tibble(Send = c(1:m)-1,Z = est_Z)
    
    pred_test_set =  test_set %>% left_join(groups, by = "Send")
    
    groups = tibble(Rec = c(1:m)-1,Z = est_Z)
    pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")
    
    baseline <- apply(non_hom_pois_est$MuA,c(1,2),mean)
    
    mu_online = rep(0,dim(pred_test_set)[1])
    for(i in 1:length(mu_online)){
      mu_online[i] = baseline[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
    }
    pred_test_set$Mean_online = mu_online
    # output_current <- pred_test_set %>% 
    #   mutate(mean_events = Mean*(end_test-end_train)) %>%
    #   mutate(diff_mean = n -mean_events) %>%
    #   mutate(diff_zero = n - 0) %>%
    #   ungroup() %>%
    #   summarise(RMSE = sqrt(mean(diff_mean^2)),
    #             RMSE_0 = sqrt(mean(diff_zero^2)))
    
    
    
    
    #-------------------------------#
    # repeat for batch estimator #
    MuA_start_batch = array(0,c(K,K,H))
    tau_start_batch = matrix(0,m,K)
    batch_nh_pois <-batch_nonhomoPois_estimator(as.matrix(train_data_batch),
                                             A_test,m,K,H,window,
                                             T=Time, dT,
                                             gravity = 0.01,
                                             MuA_start_batch,tau_start_batch,               
                                             itermax = 50,
                                             stop_eps = 0.1)
    est_Z_batch <- apply(batch_nh_pois$tau,1,which.max)
    
    groups = tibble(Send = c(1:m)-1,Z = est_Z_batch)
    pred_test_set =  pred_test_set %>% left_join(groups, by = "Send")
    
    groups = tibble(Rec = c(1:m)-1,Z = est_Z_batch)
    pred_test_set = pred_test_set %>% left_join(groups,by = "Rec")
    baseline_batch <- apply(batch_nh_pois$MuA,c(1,2),mean)
    
    mu_batch = rep(0,dim(pred_test_set)[1])
    for(i in 1:length(mu_batch)){
      mu_batch[i] = baseline_batch[pred_test_set$Z.x.x[i],pred_test_set$Z.y.y[i]]
    }
    pred_test_set$Mean_batch = mu_batch
    
    
    
    # -------- combine the results from both ---- #
    pred_test_set_events <- pred_test_set %>% 
      mutate(mean_events_online = Mean_online*(end_test-end_train),
             mean_events_batch = Mean_batch*(end_test-end_train))
    
    
    out_pred <- pred_test_set_events %>% group_by(Send) %>%
      mutate(true_out = sum(n),est_online = sum(mean_events_online),
             est_batch = sum(mean_events_batch)) %>%
      dplyr::select(Send,true_out,est_online,est_batch) %>%
      distinct(,.keep_all = T) %>%
      ungroup() %>%
      summarise(mean_resid_online = 
                  mean( (true_out-est_online)/sqrt(est_online) ),
                mean_resid_batch = 
                  mean( (true_out-est_batch)/sqrt(est_batch) )
      )
      # summarise(RMSE_online = sqrt(mean((true_out-est_online)^2)),
      #           RMSE_batch = sqrt(mean((true_out-est_batch)^2)))
    
    # print(output_current)
    # output_current$end_train <- end_train
    # output_current$sim <- sims
    # output <- output %>% bind_rows(output_current)
    out_pred$end_train <- end_train
    out_pred$sim <- sims
    output_edge_nh <- output_edge_nh %>% bind_rows(out_pred)
    
  }
  
}
beepr::beep()

# now would need to pivot_longer to compare online and batch here


output_edge_nh %>% 
  pivot_longer(-c(end_train,sim),names_to = "Model",values_to = "res") %>%
  ggplot(aes(end_train,res,color=Model)) + geom_smooth() +
  xlab("Hours history") + 
  ylab("Mean residual") +
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22)) + 
  scale_color_hue(labels=c("batch","online"))+
  #theme_bw() + 
  theme(axis.title = element_text(size=14))
