#### Real Data Analysis

## need to do online loss, and online link prediction for this here

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
##library(ppsbm)
source(here("functions/", "utils.R"))
source(here("functions/init_fcn.R"))

set.seed(100)

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid
### Read in raw data and preprocess to our format ####
# this section uses dplyr
data <- read.csv(gzfile(here("Data/sx-mathoverflow.txt.gz")))

# dim(college)

colnames(data) <- c("Data")


math_data <- data %>%
  separate(Data, c("Send", "Rec", "Time"), sep = " ") %>% 
  mutate(Time = as.numeric(Time)) %>% 
  mutate(Time = Time - min(Time)) %>% 
  mutate(Time = Time/(3600*24)) %>% 
  ## weeks instead of days
  filter(Send != Rec) %>% 
  arrange(Time)

## in the right format now

users_math <- unique(c(math_data$Send, math_data$Rec))
m_math <- length(users_math)
m_math # same as on website

math_data <- math_data %>% 
  mutate(Send = as.numeric(factor(Send,levels = users_math)) - 1) %>%
  mutate(Rec = as.numeric(factor(Rec,levels = users_math)) - 1)

## then construct A

A_math = list()


for(i in 1:nrow(math_data)){
  j = math_data$Send[i]
  k = math_data$Rec[i]
  #print(j)
  A_math[[j+1]] = j
  A_math[[k+1]] = k
}

for(i in 1:nrow(math_data)){
  j = math_data$Send[i]
  k = math_data$Rec[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A_math[[j+1]] = c(A_math[[j+1]], math_data$Rec[i])
}


A_test_math <- lapply(A_math, unique)

math_train <- math_data %>% 
  top_frac(-.85)
math_test <- math_data %>% 
  top_frac(.15)

### most events happen early on
math_train_time <- 1914
math_test_time <- 2350


math_test_set <- math_test %>% 
  group_by(Send,Rec) %>% 
  tally()


#### Hom Poisson ####
dT <- 6 # for math
K <- 3 # 2 for college, 4 for email, 3 for math
Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
Mu <- matrix(runif(K * K), K, K)
tau <- matrix(runif(m_math * K), nrow = m_math, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(0, nrow = m_math, ncol = K)

# online estimate
a <- bench::mark(
  results_pois_train <- estimate_Poisson_minimal(full_data = as.matrix(math_train),
                                         A_test_math,
                                         m_math,
                                         K,
                                         T = math_train_time,
                                         dT,
                                         # step_size = 1,
                                         B),
  iterations = 1)

# batch estimate for hom-Poisson needed here
b <- bench::mark(
  results_pois_batch <- batch_estimator_hom_Poisson(
    alltimes = as.matrix(math_train),
    A_test_math,
    m_math,
    K,
    T = math_train_time,
    itermax = 100,
    stop_eps = 0.01),
  iterations = 1
)

# simulate events
est_Z <- apply(results_pois_train$tau, 1, which.max)

groups <- tibble(Send = c(1:m_math) - 1, Z = est_Z)

pred_test_set <- math_test_set %>% 
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m_math) - 1, Z = est_Z)
pred_test_set <- pred_test_set %>% 
  left_join(groups,by = "Rec")



mu <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] <- results_pois_train$B[pred_test_set$Z.x[i], pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu


online_poss_pred <- pred_test_set %>% 
  mutate(mean_events = Mean*(math_test_time - math_train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

#### batch poisson predictions
est_Z_batch <- apply(results_pois_batch$tau, 1, which.max)

groups <- tibble(Send = c(1:m_math) - 1, Z = est_Z_batch)

pred_test_set_batch <- math_test_set %>% 
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m_math) - 1, Z = est_Z_batch)
pred_test_set_batch <- pred_test_set_batch %>% 
  left_join(groups,by = "Rec")



mu_batch <- rep(0, dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu_batch[i] <- results_pois_batch$B[pred_test_set$Z.x[i],
                                      pred_test_set$Z.y[i]]
}

pred_test_set_batch$Mean <- mu_batch


batch_poss_pred <- pred_test_set_batch %>% 
  mutate(mean_events = Mean * ( math_test_time - math_train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

pois_results <- bind_rows(online_poss_pred, batch_poss_pred) %>% 
  mutate(time = c(a$total_time, b$total_time),
         model = c("Online Poisson", "Batch Poisson"))


### Hom Hawkes ####
K <- 3 # 4 for email, 2 for college, 3 for math
dT <- 6  # 2 for email, 0.5 for college, 6 for math
Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
Mu <- matrix(runif(K * K), K, K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau <- matrix(1/K, nrow = m_math, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m_math, ncol = K)

# online estimator
## this is taking ages
a <- bench::mark(
  results_hawkes_sim <- online_estimator_eff_revised(as.matrix(math_train), 
                                                     A_test_math,
                                                     m_math,
                                                     K,
                                                     T = math_train_time,
                                                     dT, 
                                                     lam = 1,
                                                     B,
                                                     Mu,
                                                     tau,
                                                     S,
                                                     inter_T = 1),
  iterations = 1)
#### this is crashing...
# batch estimator..
b <- bench::mark(
  results_hawkes_batch <- batch_estimator(as.matrix(math_train), 
                                          A_test_math,
                                          m_math,
                                          K,
                                          T = math_train_time,
                                          dT,
                                          lam = 0.001,
                                          B,
                                          Mu,
                                          tau,
                                          itermax = 100,
                                          stop_eps = 0.01),
  iterations = 1)



est_Z <- apply(results_hawkes_sim$tau, 1, which.max)


est <- Predict_Counts_Hawkes(startT = math_train_time,
                             finalT = math_test_time,
                             A_test_math,
                             est_Z - 1,
                             results_hawkes_sim$Mu,
                             results_hawkes_sim$B,
                             m_math,
                             results_hawkes_sim$lam)

est <- est[-1,]
colnames(est) <- c("Send","Rec","Time")
est <- as_tibble(est)

hawkes_pred <- math_test_set %>% left_join(est,by = c("Send","Rec")) %>%
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


est_batch <- Predict_Counts_Hawkes(startT = math_train_time,
                                   finalT = math_test_time,
                                   A_test_math,
                                   est_Z_batch - 1,
                                   results_hawkes_batch$Mu,
                                   results_hawkes_batch$B,
                                   m_math,
                                   results_hawkes_batch$lam)

est_batch <- est_batch[-1,]
colnames(est_batch) <- c("Send", "Rec", "Time")
est_batch <- as_tibble(est_batch)

hawkes_batch_pred <- math_test_set %>%
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
  mutate(time = c(a$total_time, b$total_time),
         model = c("Online Hawkes", "Batch Hawkes"))


#### Inhom Poisson ####
window <- 1
K <- 3 # 4 for email, 2 for college, 3 for math
H <- 7
dT <- 6 # 2 for email, 0.5 for college, 6 for math
MuA_start <- array(runif(K*K*H), c(K, K, H))
tau_start <- matrix(1/K, m_math, K)
a <- bench::mark(
  non_hom_pois_est <- nonhomoPois_estimator(as.matrix(math_train),
                                            A_test_math,
                                            m_math,
                                            K,
                                            H,
                                            window,
                                            T = math_train_time,
                                            dT,
                                            gravity = 0.001,
                                            MuA_start,
                                            tau_start),
  iterations = 1)

b <- bench::mark(
  batch <- batch_nonhomoPois_estimator(as.matrix(math_train),
                                       A_test_math,
                                       m_math,
                                       K,
                                       H,
                                       window,
                                       T = math_train_time,
                                       dT,
                                       gravity = 0.001,
                                       MuA_start,
                                       tau_start,
                                       itermax = 100,
                                       stop_eps = 0.01 ),
  iterations = 1)

# taking the average of these basis functions for link prediction
# dim(non_hom_pois_est$MuA)

baseline <- apply(non_hom_pois_est$MuA, c(1, 2), mean)
est_Z <- apply(non_hom_pois_est$tau, 1, which.max)

groups <- tibble(Send = c(1:m_math) - 1, Z = est_Z)

pred_test_set <- math_test_set %>%
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m_math) - 1, Z = est_Z)
pred_test_set <- pred_test_set %>% left_join(groups, by = "Rec")



mu <- rep(0, dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] <- baseline[pred_test_set$Z.x[i], pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu

in_pois_pred <- pred_test_set %>% 
  mutate(mean_events = Mean*(math_test_time - math_train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

## batch estimates
baseline_batch <- apply(batch$MuA, c(1, 2), mean)
est_Z_batch <- apply(batch$tau, 1, which.max)

groups_batch <- tibble(Send = c(1:m_math) - 1, Z = est_Z_batch)

pred_test_set <- math_test_set %>% 
  left_join(groups_batch, by = "Send")

groups_batch <- tibble(Rec = c(1:m_math) - 1, Z = est_Z_batch)
pred_test_set <- pred_test_set %>% 
  left_join(groups, by = "Rec")



mu_batch <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu_batch[i] <- baseline[pred_test_set$Z.x[i], pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu_batch

in_pois_batch <- pred_test_set %>% 
  mutate(mean_events = Mean*(math_test_time - math_train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

in_pois_results <- bind_rows(in_pois_pred,
                             in_pois_batch) %>% 
  mutate(time = c(a$total_time, b$total_time),
         model = c("Online InPois", "Batch InPois"))


#### Inhom Hawkes ####

B_start <- matrix(runif(K * K), K, K)
a <- bench::mark(
  non_homo_Hawkes_est <- nonhomoHak_estimator_eff_revised(
    as.matrix(math_train),
    A_test_math,
    m_math,
    K,
    H,
    window,
    lam = 1,
    T = math_train_time,
    dT,
    gravity = 0.001), iterations = 1 )

## batch estimator
b <- bench::mark(
  non_homo_Hawkes_batch <- batch_nonhomoHak_estimator(
    as.matrix(math_train),
    A_test_math,
    m_math,
    K,
    H,
    window,
    T = math_train_time,
    dT,
    gravity = 0.001,
    MuA_start = MuA_start,
    tau_start = tau_start,
    lam = 0.001,
    B_start = B_start,
    itermax = 100,
    stop_eps = 0.01), iterations = 1)



# non_homo_Hawkes_est$MuA
# non_homo_Hawkes_est$Pi
# non_homo_Hawkes_est$B


baseline <- apply(non_homo_Hawkes_est$MuA, c(1, 2), mean)

est_Z <- apply(non_homo_Hawkes_est$tau, 1, which.max)
# est_Z

est <- Predict_Counts_Hawkes(startT = math_train_time,
                             finalT = math_test_time,
                             A_test_math,
                             est_Z - 1,
                             baseline,
                             non_homo_Hawkes_est$B,
                             m_math,
                             non_homo_Hawkes_est$lam)

est <- est[-1,]
colnames(est) <- c("Send","Rec","Time")
est <- as_tibble(est)

in_hawkes_pred <- math_test_set %>%
  left_join(est, by = c("Send", "Rec")) %>%
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

est_batch <- Predict_Counts_Hawkes(startT = math_train_time,
                                   finalT = math_test_time,
                                   A_test_math,
                                   est_Z_batch - 1,
                                   baseline_batch,
                                   non_homo_Hawkes_batch$B,
                                   m_math,
                                   non_homo_Hawkes_batch$lam)

est_batch <- est_batch[-1, ]
colnames(est_batch) <- c("Send","Rec","Time")
est_batch <- as_tibble(est_batch)

in_hawkes_batch <- math_test_set %>%
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
  mutate(time = c(a$total_time,
                  b$total_time),
         model = c("Online InHawkes", "Batch InHawkes"))


results <- bind_rows(pois_results,
                     hawkes_results,
                     in_pois_results,
                     in_hawkes_results) %>% 
  mutate(sim = sim_id)

results

saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("Math_Link_pred_",
                                    sim_id, ".RDS")))
