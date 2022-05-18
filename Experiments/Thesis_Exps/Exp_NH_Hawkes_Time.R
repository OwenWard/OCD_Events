#### Poisson Simulation
.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)
source(here("Experiments/", "utils.R"))
# source(here("functions/init_fcn.R"))
## Compare log likelihood, and time

no_sims <- 50

sparsity <- 0.1
Time <- 100
dT <- 1
K <- 2
H <- 7
MuA <- array(0, c(K, K, H))
MuA[, , 1] <- matrix(c(0.8, 0.2, 0.3, 0.6), 2, 2, byrow = T)
MuA[, , 2] <- matrix(c(0, 0.5, 0.2, 0.7), 2, 2, byrow = T)
MuA[, , 3] <- matrix(c(1.0, 0, 0.5, 0.6), 2, 2, byrow = T)
MuA[, , 4] <- matrix(c(0.4, 0.4, 0, 0.7), 2, 2, byrow = T)
MuA[, , 5] <- matrix(c(0.5, 0.5, 0.5, 0), 2, 2, byrow = T)
MuA[, , 6] <- matrix(c(0.4, 0, 0.4, 0), 2, 2, byrow = T)
MuA[, , 7] <- matrix(c(0, 0.9, 0, 0.9), 2, 2, byrow = T)
B <- matrix(0, K, K, byrow = TRUE)
diag(B) <- 0.5
m <- 100
Pi <- matrix(c(0.6, 0.4), 1, K)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]))
window <- 1 / 7


A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  num_edge = m * sparsity
  edge <- sample(m, num_edge) - 1
  edge <- sort(edge[edge!= i-1]) ## remove self edges
  A[[i]] <- edge
}

system.time(alltimes <- sampleBlockHak_nonhomo(Time,
                                               A, Z, MuA, B, window,
                                               lam = 1))


results <- tibble()

for(sim in 1:no_sims){
  cat("Sim:", sim, "-----\n")
  init_tau <- matrix(1/K, nrow = m, ncol = K)
  init_MuA <- array(runif(K * K * H), dim = c(K, K, H))
  init_B <- matrix(runif(K * K), nrow = K, ncol = K)
  
  
  t1 <- bench::mark(
    results.online <- nonhomoHak_estimator(alltimes,
                                           A,
                                           m, K, H, window, Time, dT,
                                           lam = 0.1,
                                           gravity = 0.001, init_B,
                                           init_MuA, init_tau
    ),
    iterations = 1)
  
  ### then compute the batch estimator
  t2 <- bench::mark(results.batch <- 
                      batch_nonhomoHak_estimator(alltimes, A, m, K, H,
                                                 window, T, dT,
                                                 lam = 0.1,
                                                 gravity = 0.001,
                                                 init_B,
                                                 init_MuA,
                                                 init_tau,
                                                 itermax = 100,
                                                 stop_eps = 0.01),
                    iterations = 1)
  ### compare the log likelihood of each
  Z_online <- apply(results.online$tau, 1, which.max) - 1
  ll_online <- get_loglik_nonhomoHak(
    alltimes, 0, Time, Z_online, results.online$MuA,
    results.online$B, results.online$Pi, A, 0.0,
    m, K, H, window
  )
  Z_batch <- apply(results.batch$tau, 1, which.max) - 1
  ll_batch <- get_loglik_nonhomoHak(
    alltimes, 0, Time, Z_batch, results.batch$MuA,
    results.online$B, results.batch$Pi, A, 0.0,
    m, K, H, window
  )
  
  ### store everything here
  diff = 1 - abs(ll_batch - ll_online)/abs(ll_batch)
  
  curr_sim <- tibble(
    sim = sim,
    t_on = t1$total_time,
    t_batch = t2$total_time,
    ll_online = ll_online,
    ll_batch = ll_batch,
    perc_diff = diff
  )
  results <- results %>% 
    bind_rows(curr_sim)
}



## then save results to something down here
saveRDS(results, file = here("Experiments",
                             "thesis_output",
                             paste0("exp_time_nh_hawkes", ".RDS")))
