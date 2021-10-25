library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("cpp_files/onlineblock.cpp")


Time <- 500
K <- 3
Mu <- matrix(c(0.6,
               0.2,
               0.3,
               0.1,
               1.0,
               0.4,
               0.5,
               0.4,
               0.8), K, K, byrow = TRUE)
#Mu <- matrix(c())
B_Pois <- matrix(0, K, K, byrow = TRUE)
m <- 100
Pi <- matrix(c(0.4, 0.3, 0.3), 1, 3)
Z <- c(rep(0, m*Pi[1]), rep(1, m*Pi[2]), rep(2, m*Pi[3]))


A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  node_list <- c(1:m)
  edge <- sample(node_list[-i], 40) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak(Time, A, Z, Mu, B_Pois, lam = 1))

nrow(alltimes)


### fit model ####
Pi <- rep(1/K,K)
B <- matrix(runif(K*K),K,K)
diag(B) <- rnorm(K,mean = 1, sd = 0.1)
tau <- matrix(runif(m*K),nrow=m,ncol=K)
tau <- tau/rowSums(tau)
S <- matrix(0,nrow = m,ncol = K)

dT <- 5
inter_T <- 5

results_online <- estimate_Poisson(full_data = alltimes,
                                   A,
                                   m,
                                   K,
                                   Time,
                                   dT,
                                   B,
                                   tau,
                                   Pi,
                                   S,
                                   inter_T,
                                   is_elbo = TRUE)


results_batch <- batch_estimator_hom_Poisson(alltimes,
                                             A,
                                             m,
                                             K,
                                             Time,
                                             B,
                                             tau,
                                             itermax = 10,
                                             stop_eps = 0.01)


### what about an even simpler problem
### only events within blocks
Time <- 500
K <- 2
Mu <- matrix(0, K, K, byrow = TRUE)
diag(Mu) <- c(1.2, 0.5)
B_Pois <- matrix(0, K, K, byrow = TRUE)
m <- 100
Pi <- matrix(c(0.5, 0.5), 1, K)
Z <- c(rep(0, m*Pi[1]), rep(1, m*Pi[2]))


A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  node_list <- c(1:m)
  # this is not good here in terms of pairs with events
  edge <- sample(node_list[-i], 40) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


system.time(alltimes <- sampleBlockHak(Time, A, Z, Mu, B_Pois, lam = 1))

nrow(alltimes)


Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
diag(B) <- rnorm(K, mean = 1, sd = 0.1)
tau <- matrix(runif(m * K), nrow = m, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m, ncol = K)

results_online <- estimate_Poisson(full_data = alltimes,
                                   A,
                                   m,
                                   K,
                                   Time,
                                   dT,
                                   B,
                                   tau,
                                   Pi,
                                   S,
                                   inter_T,
                                   is_elbo = TRUE)
results_online$B
plot(results_online$AveELBO, type = "l")


## batch estimation
results_batch <- batch_estimator_hom_Poisson(alltimes,
                                             A,
                                             m,
                                             K,
                                             Time,
                                             B,
                                             tau,
                                             itermax = 20,
                                             stop_eps = 0.01)
results_batch$tau




#### new simulation function ####

library(ppsbm)
## use this procedure for better simulation

## need piecewise const intensities
intens1 <- c(1.2)
intens2 <- c(0.75)

intens <- matrix(c(intens1, 0, 0.1, intens2), 4, 1)
# this is essentially the K*K matrix stretched out

system.time(sim1 <- generateDynppsbmConst(intens = intens, Time = 100,
                      n = 100, prop.groups = c(0.5, 0.5)))
### then need to covert this to our format


length(sim1$data$time.seq)
sim1$data$type.seq # need to convert this to pairs

# need the opposite of convert node pair

all_pairs <- listNodePairs(n = 100)
all_pairs_ind <- cbind(all_pairs, 1:nrow(all_pairs))
## this is what we need to convert back to our format, using a join
## I guess?

# this actually gives it directly too
a <- all_pairs[sim1$data$type.seq,]


as_tibble(a) %>% 
  rename(send = V1, rec = V2) %>% 
  group_by(send, rec) %>% 
  tally() %>% 
  nrow()
# so now every node pair has an event


sim_trip <- cbind(a, sim1$data$time.seq)
sim_trip


K <- 2
m <- 100
Pi <- rep(1/K, K)
B <- matrix(runif(K * K), K, K)
diag(B) <- rnorm(K, mean = 1, sd = 0.1)
tau <- matrix(runif(m * K), nrow = m, ncol = K)
tau <- tau/rowSums(tau)
S <- matrix(1/K, nrow = m, ncol = K)
### need to update the A here
tidy_sim <- as_tibble(sim_trip) %>% 
  rename(send = V1, rec = V2, time = V3)
A_new <- list()
for(i in 1:m){
  edge <- tidy_sim %>% 
    filter(send == i) %>% 
    distinct(rec) %>% 
    arrange(rec) %>% 
    pull(rec)
  A_new[[i]] <- edge - 1
}

sim_trip[, 1] <- sim_trip[, 1] - 1
sim_trip[, 2] <- sim_trip[, 2] - 1

results_online <- estimate_Poisson(full_data = sim_trip,
                                   A_new,
                                   m,
                                   K,
                                   Time,
                                   dT,
                                   B,
                                   tau,
                                   Pi,
                                   S,
                                   inter_T,
                                   is_elbo = TRUE)
results_online$B
plot(results_online$AveELBO, type = "l")


results_batch <- batch_estimator_hom_Poisson(sim_trip,
                                             A_new,
                                             m,
                                             K,
                                             Time,
                                             B,
                                             tau,
                                             itermax = 20,
                                             stop_eps = 0.01)
results_batch$tau
