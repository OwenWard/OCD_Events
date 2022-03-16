### attempting to run some simple exp to see how well we 
### can do in some simple scenarios here

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

library(kernlab)
library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/", "df_to_adj.R"))
source(here("functions", "pensky_fcns.R"))



no_sims <- 50
Time <- 100
all_results <- tibble()


for(sim in 1:no_sims){
  cat(sim, "-----\n")
  curr_sim_result <- tibble()
  
  n <- 50
  prop.groups <- c(0.5, 0.5)
  
  
  intens <- list(NULL)
  intens[[1]] <- list(intens = function(x) abs(2*sin(x/8) + 4),
                      max = 4)
  
  # curve(2*((x-7)/2), from = 0, to = 100)
  
  
  # (q,l) = (1,1)
  intens[[2]] <- list(intens = function(x) (x<=50)*0.2*(x-7)/2 +
                        (x>50)*(10*exp(-2*x/10) + 3),
                      max = 5)
  # curve( (x<=50)*0.2*(x-7)/2 +
  #          (x>50)*(10*exp(-2*x/10) + 3),
  #        from = 0, to = 100)
  
  # (q,l) = (1,2)
  intens[[3]] <- list(intens = function(x) (x<60)*2 +
                        (x>=60)*3*abs(cos(x/6) + 1),
                      max = 6)
  
  # curve( (x<60)*2 +
  #          (x>=60)*3*abs(cos(x/6) + 1), from = 0, to = 100)
  # (q,l) = (2,2)
  ## below only for directed
  intens[[4]] <- list(intens = function(x) 1 + 
                        (exp(-9*abs(x/12-1)) - .049 + 0.23*cos(x/2)),
                      max = 8)
  
  # curve(1 + 
  #         (exp(-9*abs(x/12-1)) - .049 + 0.23*cos(x/2)), from = 0, to = 100)
  
  # generate data :
  obs <- generateDynppsbm(intens,
                          Time = Time,
                          n = n,
                          prop.groups,
                          directed = TRUE)
  # latent variables (true clustering of the individuals)
  # obs$z
  # length(unique(obs$data$type.seq)) ## if n*(n-1) then all nodes have events
  
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
  ### this works
  
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
  ## this also works
  
  
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
  MuA <- array(runif(K * K * H), c(K, K, H))
  window <- 0.5
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
  
  
  ari_in_pois <- aricode::ARI(apply(results_online_inpois$tau, 1, which.max),
                              z_true)
  
  ## Fitting inhomogeneous Hawkes
  
  window <- 1/2
  K <- 2 # 4 for email, 2 for college, 3 for math
  H <- 2
  #### this was the bug#####
  dT <- 1 # 2 for email, 0.5 for college, 6 for math
  MuA_start <- array(runif(K * K * H), c(K, K, H))
  tau_start <- matrix(1/K, m, K)
  B_start <- matrix(runif(K * K), nrow = K, ncol = K)
  
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
                                                   MuA = MuA_start,
                                                   tau = tau_start)
  
  
  ari_in_haw <- aricode::ARI(z_true, apply(result_inHaw$tau, 1, which.max))
  
  ### then store all this
  curr_sim_result$ARI <- c(sc_ari, pz_ari,
                           clust_est,
                           haw_clust_est,
                           ari_in_pois,
                           ari_in_haw)
  curr_sim_result$Method <- c("Count", "PZ", "Poisson",
                              "Hawkes", "InPois", "InHawkes")
  all_results <- all_results %>% bind_rows(curr_sim_result)
}


saveRDS(all_results, file = here("Experiments",
                             "exp_results",
                             "fig_1_exp_1.RDS"))

