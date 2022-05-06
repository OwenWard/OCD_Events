#### Real Data Analysis

## need to do online loss, and online link prediction for this here

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))

set.seed(100)

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid
### Read in raw data and preprocess to our format ####
# this section uses dplyr
data <- read.csv(gzfile(here("Data/email-Eu-core-temporal.txt.gz")))

# dim(college)

colnames(data) <- c("Data")


email_data <- data %>%
  separate(Data, c("Send", "Rec", "Time"), sep = " ") %>% 
  mutate(Time = as.numeric(Time)) %>% 
  mutate(Time = Time - min(Time)) %>% 
  mutate(Time = Time/(3600*24)) %>% 
  filter(Send != Rec)

## in the right format now

users_email <- unique(c(email_data$Send, email_data$Rec))
m_email <- length(users_email)
m_email # same as on website

email_data <- email_data %>% 
  mutate(Send = as.numeric(factor(Send, levels = users_email)) - 1) %>%
  mutate(Rec = as.numeric(factor(Rec, levels = users_email)) - 1)

## then construct A

A_email <- list()


for(i in 1:nrow(email_data)){
  j = email_data$Send[i]
  k = email_data$Rec[i]
  #print(j)
  A_email[[j+1]] = j
  A_email[[k+1]] = k
}

for(i in 1:nrow(email_data)){
  j = email_data$Send[i]
  k = email_data$Rec[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A_email[[j+1]] = c(A_email[[j+1]], email_data$Rec[i])
}


A_test_email <- lapply(A_email, unique)

email_train <- email_data %>% 
  top_frac(-.85)
email_test <- email_data %>% 
  top_frac(.15)

### most events happen early on
email_train_time <- 470.3
email_test_time <- 803.9


email_test_set <- email_test %>% 
  group_by(Send,Rec) %>% 
  tally()


#### Hom Poisson ####
dT <- 2 # for math
K <- 4 # 2 for college, 4 for email, 3 for math
Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
Mu <- matrix(runif(K * K), K, K)
tau <- matrix(runif(m_email * K), nrow = m_email, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(0, nrow = m_email, ncol = K)

# online estimate
a <- system.time(
  results_pois_train <- estimate_Poisson(full_data = as.matrix(email_train),
                                         A_test_email,
                                         m_email,
                                         K,
                                         T = email_train_time,
                                         dT,
                                         B,
                                         inter_T = 5) )


b <- system.time(
  results_pois_batch <- batch_estimator_hom_Poisson(
    alltimes = as.matrix(email_train),
    A_test_email,
    m_email,
    K,
    T = email_train_time,
    itermax = 50,
    stop_eps = 0.01)
)

# simulate events
est_Z <- apply(results_pois_train$tau, 1, which.max)

groups <- tibble(Send = c(1:m_email) - 1, Z = est_Z)

pred_test_set <- email_test_set %>% 
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m_email) - 1, Z = est_Z)
pred_test_set <- pred_test_set %>% 
  left_join(groups,by = "Rec")



mu <- rep(0, dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] <- results_pois_train$B[pred_test_set$Z.x[i],
                                pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu


online_poss_pred <- pred_test_set %>% 
  mutate(mean_events = Mean*(email_test_time - email_train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

#### batch poisson predictions
est_Z_batch <- apply(results_pois_batch$tau, 1, which.max)

groups <- tibble(Send = c(1:m_email) - 1, Z = est_Z_batch)

pred_test_set_batch <- email_test_set %>% 
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m_email) - 1, Z = est_Z_batch)
pred_test_set_batch <- pred_test_set_batch %>% 
  left_join(groups,by = "Rec")



mu_batch <- rep(0, dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu_batch[i] <- results_pois_batch$B[pred_test_set$Z.x[i],
                                      pred_test_set$Z.y[i]]
}

pred_test_set_batch$Mean <- mu_batch


batch_poss_pred <- pred_test_set_batch %>% 
  mutate(mean_events = Mean*(email_test_time - email_train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

pois_results <- bind_rows(online_poss_pred, batch_poss_pred) %>% 
  mutate(time = c(as.numeric(a)[3], as.numeric(b)[3]),
         model = c("Online Poisson", "Batch Poisson"))


### Hom Hawkes ####
K <- 4 # 4 for email, 2 for college, 3 for math
dT <- 2  # 2 for email, 0.5 for college, 6 for math
Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
Mu <- matrix(runif(K * K), K, K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau <- matrix(1/K, nrow = m_email, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m_email, ncol = K)

# online estimator
## this is taking ages
a <- system.time(
  results_hawkes_sim <- online_estimator_eff_revised(as.matrix(email_train), 
                                                     A_test_email,
                                                     m_email,
                                                     K,
                                                     T = email_train_time,
                                                     dT, 
                                                     lam = 1,
                                                     B,
                                                     Mu,
                                                     tau,
                                                     inter_T = 1))
#### this is crashing...
# batch estimator..
b <- system.time(
  results_hawkes_batch <- batch_estimator(as.matrix(email_train), 
                                          A_test_email,
                                          m_email,
                                          K,
                                          T = email_train_time,
                                          dT,
                                          lam = 0.00001,
                                          B,
                                          Mu,
                                          tau,
                                          itermax =  50,
                                          stop_eps = 0.01))



est_Z <- apply(results_hawkes_sim$tau, 1, which.max)


est <- Predict_Counts_Hawkes(startT = email_train_time,
                             finalT = email_test_time,
                             A_test_email,
                             est_Z - 1,
                             results_hawkes_sim$Mu,
                             results_hawkes_sim$B,
                             m_email,
                             results_hawkes_sim$lam)

est <- est[-1,]
colnames(est) <- c("Send","Rec","Time")
est <- as_tibble(est)

hawkes_pred <- email_test_set %>% left_join(est, by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean - num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))

## batch predictions
est_Z_batch <- apply(results_hawkes_batch$tau, 1, which.max)


est_batch <- Predict_Counts_Hawkes(startT = email_train_time,
                                   finalT = email_test_time,
                                   A_test_email,
                                   est_Z_batch - 1,
                                   results_hawkes_batch$Mu,
                                   results_hawkes_batch$B,
                                   m_email,
                                   results_hawkes_batch$lam)

est_batch <- est_batch[-1,]
colnames(est_batch) <- c("Send", "Rec", "Time")
est_batch <- as_tibble(est_batch)

hawkes_batch_pred <- email_test_set %>%
  left_join(est_batch, by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))

hawkes_results <- bind_rows(hawkes_pred,
                            hawkes_batch_pred) %>% 
  mutate(time = c(as.numeric(a)[3], as.numeric(b)[3]),
         model = c("Online Hawkes", "Batch Hawkes"))


#### Inhom Poisson ####
window <- 1/7
K <- 4 # 4 for email, 2 for college, 3 for math
H <- 7
dT <- 2 # 2 for email, 0.5 for college, 6 for math
MuA_start <- array(runif(K * K * H), c(K, K, H))
tau_start <- matrix(1/K, m_email, K)
a <- system.time(
  non_hom_pois_est <- nonhomoPois_estimator(as.matrix(email_train),
                                            A_test_email,
                                            m_email,
                                            K,
                                            H,
                                            window,
                                            T = email_train_time,
                                            dT,
                                            gravity = 0.001,
                                            MuA_start,
                                            tau_start) )

b <- system.time(
  batch <- batch_nonhomoPois_estimator(as.matrix(email_train),
                                       A_test_email,
                                       m_email,
                                       K,
                                       H,
                                       window,
                                       T = email_train_time,
                                       dT,
                                       gravity = 0.001,
                                       MuA_start,
                                       tau_start,
                                       itermax = 50,
                                       stop_eps = 0.01 ))

# taking the average of these basis functions for link prediction
# dim(non_hom_pois_est$MuA)

baseline <- apply(non_hom_pois_est$MuA, c(1, 2), mean)
est_Z <- apply(non_hom_pois_est$tau, 1, which.max)

groups <- tibble(Send = c(1:m_email) - 1, Z = est_Z)

pred_test_set <- email_test_set %>% left_join(groups, by = "Send")

groups = tibble(Rec = c(1:m_email) - 1, Z = est_Z)
pred_test_set <- pred_test_set %>% left_join(groups, by = "Rec")



mu <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] <- baseline[pred_test_set$Z.x[i], pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu

in_pois_pred <- pred_test_set %>% 
  mutate(mean_events = Mean*(email_test_time - email_train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

## batch estimates
baseline_batch <- apply(batch$MuA, c(1, 2), mean)
est_Z_batch <- apply(batch$tau, 1, which.max)

groups_batch <- tibble(Send = c(1:m_email) - 1, Z = est_Z_batch)

pred_test_set <- email_test_set %>% 
  left_join(groups_batch, by = "Send")

groups_batch <- tibble(Rec = c(1:m_email) - 1, Z = est_Z_batch)
pred_test_set <- pred_test_set %>% 
  left_join(groups, by = "Rec")



mu_batch <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu_batch[i] <- baseline[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu_batch

in_pois_batch <- pred_test_set %>% 
  mutate(mean_events = Mean*(email_test_time - email_train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

in_pois_results <- bind_rows(in_pois_pred,
                             in_pois_batch) %>% 
  mutate(time = c(as.numeric(a)[3], as.numeric(b)[3]),
         model = c("Online InPois", "Batch InPois"))


#### Inhom Hawkes ####

B_start <- matrix(runif(K * K), K, K)
a <- system.time(
  non_homo_Hawkes_est <- nonhomoHak_estimator_eff_revised(
    as.matrix(email_train),
    A_test_email,
    m_email,
    K,
    H,
    window,
    lam = 1,
    T = email_train_time,
    dT,
    gravity = 0.001,
    MuA_start = MuA_start,
    tau_start = tau_start,
    B_start = B_start) )

## batch estimator
b <- system.time(
  non_homo_Hawkes_batch <- batch_nonhomoHak_estimator(
    as.matrix(email_train),
    A_test_email,
    m_email,
    K,
    H,
    window,
    T = email_train_time,
    dT,
    gravity = 0.001,
    MuA_start = MuA_start,
    tau_start = tau_start,
    lam = 0.001,
    B_start = B_start,
    itermax = 50,
    stop_eps = 0.01))


baseline <- apply(non_homo_Hawkes_est$MuA, c(1, 2), mean)

est_Z <- apply(non_homo_Hawkes_est$tau,1,which.max)
est_Z

est <- Predict_Counts_Hawkes(startT = email_train_time,
                             finalT = email_test_time,
                             A_test_email,
                             est_Z - 1,
                             baseline,
                             non_homo_Hawkes_est$B,
                             m_email,
                             non_homo_Hawkes_est$lam)

est <- est[-1,]
colnames(est) <- c("Send","Rec","Time")
est <- as_tibble(est)

in_hawkes_pred <- email_test_set %>% 
  left_join(est,by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))

## batch link prediction

baseline_batch <- apply(non_homo_Hawkes_batch$MuA, c(1, 2), mean)

est_Z_batch <- apply(non_homo_Hawkes_batch$tau, 1, which.max)

est_batch <- Predict_Counts_Hawkes(startT = email_train_time,
                                   finalT = email_test_time,
                                   A_test_email,
                                   est_Z_batch - 1,
                                   baseline_batch,
                                   non_homo_Hawkes_batch$B,
                                   m_email,
                                   non_homo_Hawkes_batch$lam)

est_batch <- est_batch[-1,]
colnames(est_batch) <- c("Send","Rec","Time")
est_batch <- as_tibble(est_batch)

in_hawkes_batch <- email_test_set %>%
  left_join(est_batch, by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))


in_hawkes_results <- bind_rows(in_hawkes_pred,
                               in_hawkes_batch) %>% 
  mutate(time = c(as.numeric(a)[3],
                  as.numeric(b)[3]),
         model = c("Online InHawkes", "Batch InHawkes"))


results <- bind_rows(pois_results,
                     hawkes_results,
                     in_pois_results,
                     in_hawkes_results) %>% 
  mutate(sim = sim_id)

results

saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("College_Link_pred_",
                                    sim_id, ".RDS")))