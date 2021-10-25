.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

nsims <- 100

results <- list()

for(sim in 1:nsims) {
  Time <- 100
  n <- 100
  intens1 <- c(1.2)
  intens2 <- c(0.75)
  true_B <- matrix(c(intens1, 0, 0.1, intens2), 
                   nrow = 2, ncol = 2, byrow = T)
  
  intens <- matrix(c(intens1, 0, 0.1, intens2), 4, 1)
  # this is essentially the K*K matrix stretched out as a single col
  
  system.time(sim1 <- generateDynppsbmConst(intens = intens,
                                            Time = Time,
                                            n = n, 
                                            prop.groups = c(0.5, 0.5)))
  
  proc_sim <- format_sims(sim_data = sim1, n = n)
  
  K <- 2
  m <- n
  Pi <- rep(1/K, K)
  B <- matrix(runif(K * K), K, K)
  # diag(B) <- rnorm(K, mean = 1, sd = 0.1)
  tau <- matrix(1, nrow = m, ncol = K)
  tau <- tau/rowSums(tau)
  S <- matrix(1/K, nrow = m, ncol = K)
  
  dT <- 1
  inter_T <- 1
  
  results_online <- estimate_Poisson(full_data = proc_sim$events,
                                     A = proc_sim$edge,
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
  
  # compute rand index
  z_true <- apply(sim1$z, 2, function(x) which(x == 1))
  z_est <- apply(results_online$tau, 1, which.max)
  clust_est <- aricode::ARI(z_true, z_est)
  
  sim_pars <- list(
    B = B, 
    est_eblo = results_online$AveELBO,
    clust = clust_est
  )
  results[[sim]] <- sim_pars
}



# #### new simulation function ####
# 
# library(ppsbm)
# ## use this procedure for better simulation
# 
# ## need piecewise const intensities
# intens1 <- c(1.2)
# intens2 <- c(0.75)
# 
# intens <- matrix(c(intens1, 0, 0.1, intens2), 4, 1)
# # this is essentially the K*K matrix stretched out
# 
# system.time(sim1 <- generateDynppsbmConst(intens = intens, Time = 100,
#                       n = 100, prop.groups = c(0.5, 0.5)))
# ### then need to covert this to our format
# 
# 
# length(sim1$data$time.seq)
# sim1$data$type.seq # need to convert this to pairs
# 
# # need the opposite of convert node pair
# 
# all_pairs <- listNodePairs(n = 100)
# all_pairs_ind <- cbind(all_pairs, 1:nrow(all_pairs))
# ## this is what we need to convert back to our format, using a join
# ## I guess?
# 
# # this actually gives it directly too
# a <- all_pairs[sim1$data$type.seq,]
# 
# 
# as_tibble(a) %>% 
#   rename(send = V1, rec = V2) %>% 
#   group_by(send, rec) %>% 
#   tally() %>% 
#   nrow()
# # so now not every node pair has an event
# 
# 
# sim_trip <- cbind(a, sim1$data$time.seq)
# sim_trip
# sim_trip[, 1] <- sim_trip[, 1] - 1
# sim_trip[, 2] <- sim_trip[, 2] - 1
# 
# K <- 2
# m <- 100
# Pi <- rep(1/K, K)
# B <- matrix(runif(K * K), K, K)
# diag(B) <- rnorm(K, mean = 1, sd = 0.1)
# tau <- matrix(runif(m * K), nrow = m, ncol = K)
# tau <- tau/rowSums(tau)
# S <- matrix(1/K, nrow = m, ncol = K)
# ### need to update the A here
# tidy_sim <- as_tibble(sim_trip) %>% 
#   rename(send = V1, rec = V2, time = V3)
# A_new <- list()
# for(i in 1:m){
#   edge <- tidy_sim %>% 
#     filter(send == i) %>% 
#     distinct(rec) %>% 
#     arrange(rec) %>% 
#     pull(rec)
#   A_new[[i]] <- edge
# }
# 
# 
# 
# Time <- 100
# dT <- 1
# inter_T <- 1
# 
# results_online <- estimate_Poisson(full_data = sim_trip,
#                                    A_new,
#                                    m,
#                                    K,
#                                    Time,
#                                    dT,
#                                    B,
#                                    tau,
#                                    Pi,
#                                    S,
#                                    inter_T,
#                                    is_elbo = TRUE)
# results_online$B
# plot(results_online$AveELBO, type = "l", ylab = "")
# 
# 
# results_batch <- batch_estimator_hom_Poisson(sim_trip,
#                                              A_new,
#                                              m,
#                                              K,
#                                              Time,
#                                              B,
#                                              tau,
#                                              itermax = 20,
#                                              stop_eps = 0.01)
# results_batch$tau
