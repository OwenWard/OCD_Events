#### Real Data Analysis

## need to do online loss, and online link prediction for this here

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
##library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/init_fcn.R"))

set.seed(100)

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid
### Read in raw data and preprocess to our format ####
# this section uses dplyr
college <- read.csv(gzfile(here("Data/CollegeMsg.txt.gz")))

# dim(college)

colnames(college) <- c("Data")


tidy_college <- college %>%
  separate(Data,c("Send","Rec","Time"), sep = " ") %>% 
  mutate(Time = as.numeric(Time)) %>% 
  mutate(Time = Time - min(Time)) %>% 
  mutate(Time = Time/(3600*24)) %>% 
  arrange(Time)

## in the right format now

users <- unique(c(tidy_college$Send, tidy_college$Rec))
m <- length(users)
m # same as on website

tidy_college <- tidy_college %>% 
  mutate(Send = as.numeric(factor(Send,levels = users))-1) %>%
  mutate(Rec = as.numeric(factor(Rec,levels = users))-1)

## then construct A

A = list()


for(i in 1:nrow(tidy_college)){
  j = tidy_college$Send[i]
  k = tidy_college$Rec[i]
  #print(j)
  A[[j+1]] = j
  A[[k+1]] = k
}

for(i in 1:nrow(tidy_college)){
  j = tidy_college$Send[i]
  k = tidy_college$Rec[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  A[[j+1]] = c(A[[j+1]], tidy_college$Rec[i])
}


A_test <- lapply(A,unique)

college_train <- tidy_college %>% 
  top_frac(-.85)
college_test <- tidy_college %>% 
  top_frac(.15)

### most events happen early on
train_time <- 76.38
test_time <- 193


test_set <- college_test %>% 
  group_by(Send,Rec) %>% 
  tally()


#### Hom Poisson ####
dT <- 1 # for math
K <- 2 # 2 for college, 4 for email, 3 for math
Pi <- rep(1/K,K)
B <- matrix(runif(K * K), K, K)
Mu <- matrix(runif(K * K), K, K)
tau <- matrix(runif(m * K), nrow = m, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(0, nrow = m, ncol = K)

# online estimate
a <- bench::mark(
  results_pois_train <- estimate_Poisson_minimal(full_data = 
                                                   as.matrix(college_train),
                                                        A_test,
                                                        m,
                                                        K,
                                                        T = train_time,
                                                        dT,
                                                        B),
  iterations = 1)

# batch estimate for hom-Poisson needed here
b <- bench::mark(
  results_pois_batch <- batch_estimator_hom_Poisson(
    alltimes = as.matrix(college_train),
    A_test,
    m,
    K,
    T = train_time,
    itermax = 100,
    stop_eps = 0.01),
  iterations = 1
)

# simulate events
est_Z <- apply(results_pois_train$tau, 1, which.max)

groups <- tibble(Send = c(1:m) - 1, Z = est_Z)

pred_test_set <- test_set %>% 
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m) - 1, Z = est_Z)
pred_test_set <- pred_test_set %>% 
  left_join(groups,by = "Rec")



mu <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] <- results_pois_train$B[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu


online_poss_pred <- pred_test_set %>% 
  mutate(mean_events = Mean*(test_time-train_time)) %>%
  mutate(diff_mean = n - mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

#### batch poisson predictions
est_Z_batch <- apply(results_pois_batch$tau, 1, which.max)

groups <- tibble(Send = c(1:m) - 1, Z = est_Z_batch)

pred_test_set_batch <- test_set %>% 
  left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m) - 1, Z = est_Z_batch)
pred_test_set_batch <- pred_test_set_batch %>% 
  left_join(groups,by = "Rec")



mu_batch <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu_batch[i] <- results_pois_batch$B[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set_batch$Mean <- mu_batch


batch_poss_pred <- pred_test_set_batch %>% 
  mutate(mean_events = Mean*(test_time-train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

pois_results <- bind_rows(online_poss_pred, batch_poss_pred) %>% 
  mutate(time = c(a$total_time, b$total_time),
         memory = c(a$mem_alloc, b$mem_alloc),
         model = c("Online Poisson", "Batch Poisson"))


### Hom Hawkes ####
K <- 2 # 4 for email, 2 for college, 3 for math
dT <- 1  # 2 for email, 0.5 for college, 6 for math
Pi <- rep(1/K,K)
B <- matrix(runif(K*K), K, K)
Mu <- matrix(runif(K*K), K, K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau <- matrix(runif(m * K), nrow = m, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m, ncol = K)

# online estimator
a <- bench::mark(
  results_hawkes_sim <- online_estimator_eff_revised(as.matrix(college_train), 
                                                     A_test,
                                                     m,
                                                     K,
                                                     T = train_time,
                                                     dT, 
                                                     lam = 1,
                                                     B,
                                                     Mu,
                                                     tau,
                                                     S,
                                                     inter_T = 1),
  iterations = 1)
#### here!
# batch estimator..
b <- bench::mark(
  results_hawkes_batch <- batch_estimator(as.matrix(college_train), 
                                          A_test,
                                          m,
                                          K,
                                          T = train_time,
                                          dT,
                                          lam = 1,
                                          B,
                                          Mu,
                                          tau,
                                          itermax = 100,
                                          stop_eps = 0.01),
  iterations = 1)



est_Z <- apply(results_hawkes_sim$tau, 1, which.max)


est <- Predict_Counts_Hawkes(startT = train_time,
                            finalT = test_time,A_test,
                            est_Z - 1,
                            results_hawkes_sim$Mu,
                            results_hawkes_sim$B,m,
                            results_hawkes_sim$lam)

est <- est[-1,]
colnames(est) <- c("Send","Rec","Time")
est <- as_tibble(est)

hawkes_pred <- test_set %>% left_join(est,by = c("Send","Rec")) %>%
  rename(Pred_Mean = Time,
         num_test = n) %>%
  mutate(Pred_Mean = replace_na(Pred_Mean,0)) %>%
  mutate(diff_mean = Pred_Mean-num_test,
         diff_zer = 0 - num_test) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zer^2)))

## batch predictions
est_Z_batch <- apply(results_hawkes_batch$tau, 1, which.max)


est_batch <- Predict_Counts_Hawkes(startT = train_time,
                             finalT = test_time,A_test,
                             est_Z_batch - 1,
                             results_hawkes_batch$Mu,
                             results_hawkes_batch$B,
                             m,
                             results_hawkes_batch$lam)

est_batch <- est_batch[-1,]
colnames(est_batch) <- c("Send", "Rec", "Time")
est_batch <- as_tibble(est_batch)

hawkes_batch_pred <- test_set %>%
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
         memory = c(a$mem_alloc, b$mem_alloc),
         model = c("Online Hawkes", "Batch Hawkes"))


#### Inhom Poisson ####
window <- 1
K <- 2 # 4 for email, 2 for college, 3 for math
H <- 7
dT <- 1 # 2 for email, 0.5 for college, 6 for math
MuA_start <- array(0, c(K, K, H))
tau_start <- matrix(1/K, m, K)
a <- bench::mark(
  non_hom_pois_est <- nonhomoPois_estimator(as.matrix(college_train),
                                            A_test,
                                            m,
                                            K,
                                            H,
                                            window,
                                            T = train_time,
                                            dT,
                                            gravity = 0.001,
                                            MuA_start,
                                            tau_start),
  iterations = 1)

b <- bench::mark(
  batch <- batch_nonhomoPois_estimator(as.matrix(college_train),
                                       A_test,
                                       m,
                                       K,
                                       H,
                                       window,
                                       T = train_time,
                                       dT,
                                       gravity = 0.001,
                                       MuA_start,
                                       tau_start,
                                       itermax = 100,
                                       stop_eps = 0.002 ),
  iterations = 1)

# taking the average of these basis functions for link prediction
# dim(non_hom_pois_est$MuA)

baseline <- apply(non_hom_pois_est$MuA, c(1, 2), mean)
est_Z <- apply(non_hom_pois_est$tau, 1, which.max)

groups <- tibble(Send = c(1:m) - 1, Z = est_Z)

pred_test_set <- test_set %>% left_join(groups, by = "Send")

groups <- tibble(Rec = c(1:m) - 1, Z = est_Z)
pred_test_set <- pred_test_set %>% left_join(groups,by = "Rec")


mu <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu[i] <- baseline[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu

in_pois_pred <- pred_test_set %>% 
  mutate(mean_events = Mean*(test_time-train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

## batch estimates
baseline_batch <- apply(batch$MuA, c(1, 2), mean)
est_Z_batch <- apply(batch$tau, 1, which.max)

groups_batch <- tibble(Send = c(1:m) - 1, Z = est_Z_batch)

pred_test_set <- test_set %>% 
  left_join(groups_batch, by = "Send")

groups_batch <- tibble(Rec = c(1:m) - 1, Z = est_Z_batch)
pred_test_set <- pred_test_set %>% 
  left_join(groups_batch, by = "Rec")



mu_batch <- rep(0,dim(pred_test_set)[1])
for(i in 1:length(mu)){
  mu_batch[i] <- baseline[pred_test_set$Z.x[i],pred_test_set$Z.y[i]]
}

pred_test_set$Mean <- mu_batch

in_pois_batch <- pred_test_set %>% 
  mutate(mean_events = Mean*(test_time-train_time)) %>%
  mutate(diff_mean = n -mean_events) %>%
  mutate(diff_zero = n - 0) %>%
  ungroup() %>%
  summarise(RMSE = sqrt(mean(diff_mean^2)),
            RMSE_0 = sqrt(mean(diff_zero^2)))

in_pois_results <- bind_rows(in_pois_pred,
                             in_pois_batch) %>% 
  mutate(time = c(a$total_time, b$total_time),
         memory = c(a$mem_alloc, b$mem_alloc),
         model = c("Online InPois", "Batch InPois"))


#### Inhom Hawkes ####

B_start <- matrix(0, K, K)
a <- bench::mark(
  non_homo_Hawkes_est <- nonhomoHak_estimator_eff_revised(
    as.matrix(college_train),
    A_test,
    m,
    K,
    H,
    window,
    lam = 1,
    T = train_time,
    dT,
    gravity = 0.001), iterations = 1 )

## batch estimator
b <- bench::mark(
  non_homo_Hawkes_batch <- batch_nonhomoHak_estimator(
    as.matrix(college_train),
    A_test,
    m,
    K,
    H,
    window,
    T = train_time,
    dT,
    gravity = 0.001,
    MuA_start = MuA_start,
    tau_start = tau_start,
    lam = 1,
    B_start = B_start,
    itermax = 100,
    stop_eps = 0.01),
  iterations = 1)



# non_homo_Hawkes_est$MuA
# non_homo_Hawkes_est$Pi
# non_homo_Hawkes_est$B


baseline <- apply(non_homo_Hawkes_est$MuA, c(1, 2), mean)

est_Z <- apply(non_homo_Hawkes_est$tau,1,which.max)
est_Z

est <- Predict_Counts_Hawkes(startT = train_time,
                             finalT = test_time,
                             A_test,
                             est_Z - 1,
                             baseline,
                             non_homo_Hawkes_est$B,
                             m,
                             non_homo_Hawkes_est$lam)

est <- est[-1,]
colnames(est) <- c("Send","Rec","Time")
est <- as_tibble(est)

in_hawkes_pred <- test_set %>% left_join(est,by = c("Send","Rec")) %>%
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

est_batch <- Predict_Counts_Hawkes(startT = train_time,
                             finalT = test_time,
                             A_test,
                             est_Z_batch - 1,
                             baseline_batch,
                             non_homo_Hawkes_batch$B,
                             m,
                             non_homo_Hawkes_batch$lam)

est_batch <- est_batch[-1,]
colnames(est_batch) <- c("Send","Rec","Time")
est_batch <- as_tibble(est_batch)

in_hawkes_batch <- test_set %>%
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
         memory = c(a$mem_alloc, b$mem_alloc),
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
