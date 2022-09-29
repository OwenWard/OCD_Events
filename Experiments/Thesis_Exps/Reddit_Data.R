#### Reddit Subreddit Data

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
data <- read.csv(here("Data/soc-redditHyperlinks-body.tsv"), sep = "\t")

data[1, ]
dim(data)
data[1:5, 1:4]

## keep the data
event_data <- data[ , 1:5]
event_data

### process this and keep the subreddit names somewhere

reddit_data <- event_data %>% 
  select("SOURCE_SUBREDDIT", "TARGET_SUBREDDIT", "TIMESTAMP", "LINK_SENTIMENT")

## convert timestamp to time since min time here
class(reddit_data$TIMESTAMP)

min_date <- min(reddit_data$TIMESTAMP)


time_sec <- as.POSIXct(reddit_data$TIMESTAMP) - (as.POSIXct(min_date) - 1)

reddit_data$time_days <- as.numeric(time_sec/(60 * 60 * 24))

reddit_data


pos_data <- reddit_data %>% filter(LINK_SENTIMENT == 1) %>% 
  filter(SOURCE_SUBREDDIT != TARGET_SUBREDDIT)


#### then do whatever and cluster this
users_reddit <- unique(c(pos_data$SOURCE_SUBREDDIT, 
                         pos_data$TARGET_SUBREDDIT))
m_reddit <- length(users_reddit)
m_reddit # same as on website

reddit_id <- tibble(ind = 1:m_reddit, name = users_reddit)

pos_data <- pos_data %>% 
  mutate(Send = as.numeric(factor(SOURCE_SUBREDDIT, 
                                  levels = users_reddit)) - 1) %>%
  mutate(Rec = as.numeric(factor(TARGET_SUBREDDIT,
                                 levels = users_reddit)) - 1) %>% 
  filter(Send != Rec)


#### 
pos_data <- pos_data %>% select(Send, Rec, time_days) %>% arrange(time_days)


### then fit to this data here
A_reddit <- list()

for(i in 1:nrow(pos_data)){
  j = pos_data$Send[i]
  k = pos_data$Rec[i]
  #print(j)
  A_reddit[[j+1]] = j
  A_reddit[[k+1]] = k
}

for(i in 1:nrow(pos_data)){
  j = pos_data$Send[i]
  k = pos_data$Rec[i]
  #print(j)
  #A[[j+1]] = c(A[[j+1]],j)
  #A[[k+1]] = c(A[[k+1]],k)
  print(k)
  A_reddit[[j+1]] = c(A_reddit[[j+1]], pos_data$Rec[i])
}


A_test_reddit <- lapply(A_reddit, unique)


dT <- 10 # for math
K <- 5 # 2 for college, 4 for email, 3 for math
# Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
# Mu <- matrix(runif(K * K), K, K)
# tau <- matrix(runif(m_math * K), nrow = m_math, ncol = K)
# tau <- tau/rowSums(tau)
# S <- matrix(0, nrow = m_math, ncol = K)

Time <- 1216
# online estimate
results_pois_train <- estimate_Poisson_minimal(full_data = as.matrix(pos_data),
                                       A_test_reddit,
                                       m_reddit,
                                       K,
                                       T = Time,
                                       dT = dT,
                                       B)


results_pois_train$Pi
results_pois_train$B
### issue here with an entry in A_test_Reddit which is empty

results_pois_batch <- batch_estimator_hom_Poisson(
  alltimes = as.matrix(pos_data),
  A_test_reddit,
  m_reddit,
  K,
  T = Time,
  itermax = 100,
  stop_eps = 0.002)

results_pois_batch$Pi
results_pois_batch$B


###### Hawkes
K <- 2 # 4 for email, 2 for college, 3 for math
dT <- 25  # 2 for email, 0.5 for college, 6 for math
Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
Mu <- matrix(runif(K * K), K, K)
#diag(B) = rnorm(3,mean = 1, sd = 0.1)
tau <- matrix(1/K, nrow = m_math, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m_math, ncol = K)

# online estimator
## this is taking ages

results_hawkes_sim <- online_estimator_eff_revised(as.matrix(pos_data), 
                                                   A_test_reddit,
                                                   m_reddit,
                                                   K,
                                                   T = Time,
                                                   dT, 
                                                   lam = 1,
                                                   B,
                                                   Mu,
                                                   tau,
                                                   inter_T = 1)

results_hawkes_sim$Pi

### hawkes just collapses to one cluster


### Inhom Poisson maybe 

window <- 1/7
K <- 4 # 4 for email, 2 for college, 3 for math
H <- 7
dT <- 1 # 2 for email, 0.5 for college, 6 for math
MuA_start <- array(runif(K*K*H), c(K, K, H))
tau_start <- matrix(1/K, m_math, K)
non_hom_pois_est <- nonhomoPois_estimator(as.matrix(pos_data),
                                            A_test_reddit,
                                            m_reddit,
                                            K,
                                            H,
                                            window,
                                            T = Time,
                                            dT,
                                            gravity = 0.001,
                                            MuA_start,
                                            tau_start)

non_hom_pois_est$Pi
