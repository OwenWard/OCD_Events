#### Experiment to generate Figure 1 in the paper 
## Illustrating the importance of incorporating the temporal
## Dynamics if present when the goal is community detection


library(here)
library(kernlab)
library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/", "df_to_adj.R"))
source(here("functions", "pensky_fcns.R"))


Time <- 100
# no_sims <- 50
dT <- 1
inter_T <- 1
K <- 2

sparsity <- 0.25 # prop of edges which can have events


### how to generate the data?

#### Simulate from PPSBM ####


n <- 50
prop.groups <- c(0.05, 0.95)
# 3 different intensity functions :
## these have no events after a short time window
# intens <- list(NULL)
# intens[[1]] <- list(intens = function(x) 100*x*exp(-8*x),
#                     max = 5)
# # (q,l) = (1,1)
# intens[[2]] <- list(intens = function(x) exp(3*x)*(sin(6*pi*x-pi/2)+1)/2,
#                     max = 13)
# # (q,l) = (1,2)
# intens[[3]] <- list(intens = function(x) 8.1*(exp(-6*abs(x-1/2))-.049),
#                     max = 8)
# (q,l) = (2,2)


intens <- list(NULL)
intens[[1]] <- list(intens = function(x) abs(2*sin(x/8) + 4)/12,
                    max = 0.5)
# (q,l) = (1,1)
intens[[2]] <- list(intens = function(x) (1 + cos(x/2))/20,
                    max = 0.01)
# (q,l) = (1,2)
intens[[3]] <- list(intens = function(x) (exp(-6*abs(x/4-1/2))-.049)/20,
                    max = 0.005)
# (q,l) = (2,2)
## below only for directed
intens[[4]] <- list(intens = function(x) 1 + (exp(-6*abs(x/4-20))-.9),
                    max = 1)


# generate data :
obs <- generateDynppsbm(intens,
                        Time = Time,
                        n = n,
                        prop.groups,
                        directed = TRUE)
# latent variables (true clustering of the individuals)
obs$z
length(unique(obs$data$type.seq)) ## if n*(n-1) then all nodes have events

### then convert this to format we can use, also
### construct aggregate adj matrices and use Pensky


### Fit SBM to count matrix

proc_sim <- format_sims(sim_data = obs, n = n, directed = TRUE)
### this doesn't seem to be working correctly for undirected

A_test <- proc_sim$edge
events <- proc_sim$events

### Fit SBM to count matrix
df <- as_tibble(proc_sim$events) %>% 
  rename(Send = V1,
         Rec = V2,
         Time = V3) %>% 
  group_by(Send, Rec) %>% 
  count()

A_mat <- summ_to_adj_mat(df, n = n)

### then do specc on this

sc <- specc(A_mat, centers = 2)
sc_est <- sc@.Data
z_true <- apply(obs$z, 2, which.max)
(sc_ari <- aricode::ARI(z_true, sc_est))
### good, this doesn't work

### fit Pensky Zhang to this data
A_mats <- event_to_mat_seq(proc_sim$events,
                           Total_time = Time, window_size = 1, n = n)
## convert this to binary
A_mats[A_mats > 0] <- 1

l <- 4
m_pz <- 4
m0 <- 1
l0 <- 2

A_pz <- pz_estimator_3(A = A_mats, time = Time,
                       l0 = l0, m0 = m0,
                       m = m_pz, r = (Time - 1))
A_pz
est_labels <- specc(A_pz, centers = 2)
(pz_ari <- aricode::ARI(z_true, est_labels@.Data))
## same poor performance


#### Fit our method to this data
K <- 2
m <- n
B <- matrix(runif(K * K), K, K)
dT <- 1
inter_T <- 1
# capture output to not print out
results_online <- estimate_Poisson(full_data = proc_sim$events,
                                   A = proc_sim$edge,
                                   m,
                                   K,
                                   Time,
                                   dT,
                                   B,
                                   inter_T,
                                   is_elbo = TRUE)

# z_true <- apply(obs$z, 2, which.max)
z_est <- apply(results_online$tau, 1, which.max)
(clust_est <- aricode::ARI(z_true, z_est))

### makes sense Poisson not that good here.


## Block Hawkes
Mu <- matrix(runif(K * K), K, K)
B <- matrix(runif(K * K), K, K)
tau <- matrix(1/K, nrow = m, ncol = K)
results_online_hawkes <- online_estimator_eff_revised(alltimes = proc_sim$events, 
                                               A = proc_sim$edge,
                                               m,
                                               K,
                                               Time,
                                               dT = dT,
                                               lam = 1,
                                               B, Mu, tau,
                                               inter_T, 
                                               is_elbo = FALSE)

z_est_haw <- apply(results_online_hawkes$tau, 1, which.max)
(haw_clust_est <- aricode::ARI(z_true, z_est_haw))

## Fitting Inhomogeneous Poisson
H <- 2
MuA <- array(0.5, c(K, K, H))
window <- 0.25
tau <- matrix(1/K, nrow = m, ncol = K)
system.time(results_online_inpois <- nonhomoPois_estimator(alltimes = events,
                                                    A = proc_sim$edge,
                                                    m,
                                                    K,
                                                    H,
                                                    window,
                                                    Time,
                                                    dT,
                                                    gravity = 0.001,
                                                    MuA, tau))

### this only works sometimes..., because window not defined


aricode::ARI(apply(results_online_inpois$tau, 1, which.max), z_true)

## Fitting inhomogeneous Hawkes

window <- 1/4
K <- 2 # 4 for email, 2 for college, 3 for math
H <- 4
dT <- 1 # 2 for email, 0.5 for college, 6 for math
MuA_start <- array(0.25, c(K, K, H))
tau_start <- matrix(1/K, m, K)
B_start <- matrix(0.5, nrow = K, ncol = K)

result_inHaw <- nonhomoHak_estimator_eff_revised(alltimes = proc_sim$events,
                         A = proc_sim$edge,
                         m,
                         K,
                         H,
                         window,
                         T = Time,
                         dT,
                         lam = 0.1,
                         gravity = 0.0,
                         B_start = B_start,
                         ### this not defined!
                         MuA = MuA_start,
                         tau = tau_start)
## this still crashing for this type of data but not for below...

apply(result_inHaw$tau, 1, which.max)
aricode::ARI(apply(result_inHaw$tau, 1, which.max), z_true)
### so in this setting we can actually all recover them, which is something


### aim is to have down to here at least working correctly ####

### simulate from Inhomogeneous Poisson ####

Time <- 100
K <- 2
H <- 2
MuA <- array(0, c(K, K, H))
MuA[,,1] <- matrix(c(0.8, 0.2, 0.6, 0.4), 2, 2)
MuA[,,2] <- matrix(c(0.4, 0.7, 0.2, 0.7), 2, 2)
B <- matrix(0, K, K, byrow = TRUE)
m <- 100
Pi <- matrix(c(0.6, 0.4), 1, K)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]))
window <- 1.5
sparsity <- 0.05

## add in some sparsity here

A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  num_edge = m * sparsity
  edge <- sample(m, num_edge) - 1
  edge <- sort(edge[edge!= i-1])
  # edge <- sort(edge)
  A[[i]] <- edge
}



system.time(alltimes <- sampleBlockHak_nonhomo(Time,
                                               A,
                                               Z,
                                               MuA,
                                               B, window, lam = 1))

dim(alltimes)


### lets fit our model to this at least
MuA_start <- array(1/K, c(K, K, H))
tau_start <- matrix(1/K, nrow = m, ncol = K)
S_start <- matrix(1/K, nrow = m, ncol = K)

est <- nonhomoPois_estimator(alltimes, A, m, K, H, window,
                             Time, dT = 1, gravity = 0.01, MuA_start,
                             tau_start, is_elbo = FALSE)

z_est <- apply(est$tau, 1, which.max)
aricode::ARI(z_est, Z)

## so can recover this one well

## what about when use count matrix?
df <- as_tibble(alltimes) %>% 
  rename(Send = V1,
         Rec = V2,
         Time = V3) %>% 
  group_by(Send, Rec) %>% 
  count()

A_mat <- summ_to_adj_mat(df, n = m)

### then do specc on this

sc <- specc(A_mat, centers = 2)
sc_est <- sc@.Data
(sc_ari <- aricode::ARI(Z, sc_est))
## so this works poorly here


## what about something over time, ala the Pensky method?
## convert this into a sequence of A matrices, one for each day?

A_mats <- event_to_mat_seq(alltimes, Total_time = Time, window_size = 1, n = m)
## convert this to binary
A_mats[A_mats > 0] <- 1

l <- 4
m_pz <- 4
m0 <- 1
l0 <- 2

A_pz <- pz_estimator_3(A = A_mats, time = Time,
                       l0 = l0, m0 = m0,
                       m = m_pz, r = (Time - 1))
A_pz
est_labels <- specc(A_pz, centers = 2)
(pz_ari <- aricode::ARI(Z, est_labels@.Data))


### try fit inhom Hawkes to this, see what happens

B_start <- matrix(0.5, nrow = K, ncol = K)
dT <- 1
system.time(results.online <- nonhomoHak_estimator_eff_revised(alltimes,
                                                   A,
                                                   m,
                                                   K,
                                                   H,
                                                   window,
                                                   T = Time,
                                                   dT,
                                                   lam = 0.1, 
                                                   gravity = 0.001,
                                                   B_start,
                                                   MuA_start,
                                                   tau_start))

aricode::ARI(Z, apply(results.online$tau, 1, which.max))

### this seems to be working now at least

hist(alltimes[,3])

### so we can recover this, but what sort of dynamics is it really?


#### need some real data, where we have some known community structure
#### or generate fake community structure for some real 
#### event data, or something like that

### ppsbm data examples
hist(generated_Q3$data$time.seq)

### this shows some temporal dynamics

###
n <- 20
Time <- 100
proc_sim <- format_sims(sim_data = generated_Q3_n20, n = n, directed = FALSE)
### this doesn't seem to be working correctly for undirected
z_true <- apply(generated_Q3_n20$z, 2, which.max)

A_test <- proc_sim$edge
events <- proc_sim$events
## scale back to [0,100] instead of [0,1]
events[,3] <- events[,3]*100

### Fit SBM to count matrix
df <- as_tibble(events) %>% 
  rename(Send = V1,
         Rec = V2,
         Time = V3) %>% 
  group_by(Send, Rec) %>% 
  count()

A_mat <- summ_to_adj_mat(df, n = n)

### then do specc on this

sc <- specc(A_mat, centers = 3)
sc_est <- sc@.Data
# z_true <- apply(obs$z, 2, which.max)
(sc_ari <- aricode::ARI(z_true, sc_est))

### then Pensky and Zhang
A_mats <- event_to_mat_seq(events,
                           Total_time = Time, window_size = 1, n = n)
## convert this to binary
A_mats[A_mats > 0] <- 1

l <- 4
m_pz <- 4
m0 <- 1
l0 <- 2

### need to check this to see what's up
A_pz <- pz_estimator_3(A = A_mats, time = Time,
                       l0 = l0, m0 = m0,
                       m = m_pz, r = (Time - 1))
A_pz
est_labels <- specc(A_pz, centers = 3)
(pz_ari <- aricode::ARI(z_true, est_labels@.Data))


### our methods
K <- 3
m <- n
B <- matrix(runif(K * K), K, K)
dT <- 0.25
inter_T <- 1
# capture output to not print out
results_online <- estimate_Poisson(full_data = events,
                                   A = proc_sim$edge,
                                   m,
                                   K,
                                   Time,
                                   dT,
                                   B,
                                   inter_T,
                                   is_elbo = TRUE)

# z_true <- apply(obs$z, 2, which.max)
z_est <- apply(results_online$tau, 1, which.max)
(clust_est <- aricode::ARI(z_true, z_est))


### Hawkes
Mu <- matrix(runif(K * K), K, K)
B <- matrix(runif(K * K), K, K)
tau <- matrix(1/K, nrow = m, ncol = K)
results_online_hawkes <- online_estimator_eff_revised(alltimes = events, 
                                                      A = proc_sim$edge,
                                                      m,
                                                      K,
                                                      Time,
                                                      dT = dT,
                                                      lam = 1,
                                                      B, Mu, tau,
                                                      inter_T, 
                                                      is_elbo = FALSE)

z_est_haw <- apply(results_online_hawkes$tau, 1, which.max)
(haw_clust_est <- aricode::ARI(z_true, z_est_haw))


### Inhom Poisson
H <- 2
dT <- 0.25
MuA <- array(0.5, c(K, K, H))
window <- 0.25
tau <- matrix(1/K, nrow = m, ncol = K)
system.time(results_online_inpois <- nonhomoPois_estimator(alltimes = events,
                                                           A = proc_sim$edge,
                                                           m,
                                                           K,
                                                           H,
                                                           window,
                                                           Time,
                                                           dT,
                                                           gravity = 0.001,
                                                           MuA, tau))

### this only works sometimes..., because window not defined


aricode::ARI(apply(results_online_inpois$tau, 1, which.max), z_true)


### Inhom Hawkes
window <- 1/8
K <- 3 # 4 for email, 2 for college, 3 for math
H <- 4
dT <- 0.25 # 2 for email, 0.5 for college, 6 for math
MuA_start <- array(runif(K*K*H), c(K, K, H))
tau_start <- matrix(1/K, m, K)
B_start <- matrix(runif(K*K), nrow = K, ncol = K)

result_inHaw <- nonhomoHak_estimator_eff_revised(alltimes = events,
                                                 A = proc_sim$edge,
                                                 m,
                                                 K,
                                                 H,
                                                 window,
                                                 T = Time,
                                                 dT,
                                                 lam = 0.1,
                                                 gravity = 0.01,
                                                 B_start = B_start,
                                                 MuA = MuA_start,
                                                 tau_start = tau_start)


z_est_inH <- apply(result_inHaw$tau, 1, which.max)
aricode::ARI(z_est_inH, z_true)

#### what about the smaller example ####
#### no better either


intens1 <- c(1,3,8)/4
intens2 <- c(2,3,6)/6
intens <- matrix(c(intens1,intens2,intens1,intens2),4,3)
Time <- 100
n <- 100
prop.groups <- c(0.5, 0.5)
dynppsbm <- generateDynppsbmConst(intens,
                                  Time,
                                  n,
                                  prop.groups,
                                  directed = TRUE)

###
hist(dynppsbm$data$time.seq)
proc_sim <- format_sims(sim_data = dynppsbm, n = n,
                        directed = TRUE)

z_true <- apply(dynppsbm$z, 2, which.max)

A_test <- proc_sim$edge
events <- proc_sim$events


### Fit SBM to count matrix
df <- as_tibble(events) %>% 
  rename(Send = V1,
         Rec = V2,
         Time = V3) %>% 
  group_by(Send, Rec) %>% 
  count()

A_mat <- summ_to_adj_mat(df, n = n)

sc <- specc(A_mat, centers = 2)
sc_est <- sc@.Data
# z_true <- apply(obs$z, 2, which.max)
(sc_ari <- aricode::ARI(z_true, sc_est))

### then Pensky and Zhang
A_mats <- event_to_mat_seq(events,
                           Total_time = Time, window_size = 1, n = n)
## convert this to binary
A_mats[A_mats > 0] <- 1

l <- 4
m_pz <- 4
m0 <- 1
l0 <- 2

### need to check this to see what's up
A_pz <- pz_estimator_3(A = A_mats, time = Time,
                       l0 = l0, m0 = m0,
                       m = m_pz, r = (Time - 1))
A_pz
est_labels <- specc(A_pz, centers = 3)
(pz_ari <- aricode::ARI(z_true, est_labels@.Data))


### our methods
K <- 2
m <- n
B <- matrix(runif(K * K), K, K)
dT <- 0.25
inter_T <- 1
# capture output to not print out
results_online <- estimate_Poisson(full_data = events,
                                   A = proc_sim$edge,
                                   m,
                                   K,
                                   Time,
                                   dT,
                                   B,
                                   inter_T,
                                   is_elbo = TRUE)

# z_true <- apply(obs$z, 2, which.max)
z_est <- apply(results_online$tau, 1, which.max)
(clust_est <- aricode::ARI(z_true, z_est))


H <- 2
dT <- 0.25
MuA <- array(0.5, c(K, K, H))
window <- 0.25
tau <- matrix(1/K, nrow = m, ncol = K)
system.time(results_online_inpois <- nonhomoPois_estimator(alltimes = events,
                                                           A = proc_sim$edge,
                                                           m,
                                                           K,
                                                           H,
                                                           window,
                                                           Time,
                                                           dT,
                                                           gravity = 0.001,
                                                           MuA, tau))

z_est_inP <- apply(results_online_inpois$tau, 1, which.max)
(clust_est_Pois <- aricode::ARI(z_true, z_est_inP))

### ok, nice! so Pensky and Zheng is a definite improvement 
### over just the counts, but maybe still not good enough
### this is definitely a good example, just need to generate it
### and simulate it


