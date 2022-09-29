##### Online Loss for Real Data


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
  mutate(Time = Time/(3600*24))

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
# train_time <- 76.38
Time <- 193



#### Hom Poisson ####

## fit batch first
K <- 10
results_pois_batch <- batch_estimator_hom_Poisson(
  alltimes = as.matrix(tidy_college),
  A_test,
  m,
  K,
  T = Time,
  itermax = 100,
  stop_eps = 0.01)

### compute best batch loss
z_batch <- apply(results_pois_batch$tau, 1, which.max)

batch_b <- results_pois_batch$B

batch_loss <- batch_loss(full_data = as.matrix(tidy_college),
                       A = A_test,
                       m,
                       K,
                       Time,
                       dT = 1, 
                       batch_B = batch_b, true_z = z_batch)
batch_average <- mean(batch_loss$Batch_loss/1:Time)


### then online loss
inter_T <- 1
dT <- 0.1
# capture output to not print out
### use initialization scheme here first also
colnames(tidy_college) <- c("V1", "V2", "V3")
result <- sparse_poisson(alltimes = tidy_college,
                        K = K,
                        n0 = 10,
                        m = m,
                        m0 = 500)

Mu_est <- result$est_B
init_tau <- matrix(0, nrow = m, ncol = K)
for(i in seq_along(result$est_clust)){
  init_tau[i, result$est_clust[i]] <- 1
}

results_online_init <- estimate_Poisson_init(full_data = 
                                               as.matrix(result$rest_events),
                                             A = A_test,
                                             m,
                                             K,
                                             Time,
                                             dT = dT,
                                             B = Mu_est,
                                             inter_T,
                                             init_tau,
                                             start = result$cut_off,
                                             is_elbo = FALSE)


B_ests <- results_online_init$inter_B
z_true <- rep(1, m)
true_B <- matrix(runif(K*K), nrow = K, ncol = K)
## need to modify the regret function also
init_time <- Time - result$cut_off
out <- compute_regret(full_data = 
                        as.matrix(result$rest_events),
                      A = A_test, 
                      m,
                      K,
                      init_time,
                      dT,
                      true_z = z_true,
                      B_ests = B_ests,
                      tau_ests = results_online_init$early_tau,
                      true_B = true_B)

## only empirical parts of this will make sense
card_A <- as_tibble(tidy_college) %>% 
  select(V1, V2) %>% 
  distinct() %>% 
  nrow()
## is this regret function correct?
est_loss <- -out$EstLLH/card_A
best_loss <- -out$TrueLLH/card_A
regret <- cumsum(est_loss) - cumsum(best_loss)

time_seq <- seq(from = Time - init_time + dT, to = Time, by = dT)
M_vec <- length(out$Online_Loss_Est) + result$cut_off 

### compute batch estimates and online loss estimates here
onl_loss <- tibble(cumsum(out$Online_Loss_True)/((result$cut_off+1):M_vec),
                   cumsum(out$Online_Loss_Est)/((result$cut_off+1):M_vec),
                   dT = time_seq) 
colnames(onl_loss) <- c("True_Z", "Est_Z", "dT")

onl_loss %>% ggplot(aes(dT, Est_Z)) + geom_line() +
  labs(x = "Time", y = "Loss") +
  geom_hline(yintercept = batch_average, colour = "red")

##### then would need to repeat this for Hawkes, etc