### attempting to run some simple exp to see how well we 
### can do in some simple scenarios here

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

library(kernlab)
library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/", "df_to_adj.R"))
source(here("functions", "pensky_fcns.R"))

print(here("Experiments", "exp_results", "fig_1_exp_1.RDS"))

no_sims <- 50
Time <- 100
all_results <- tibble()


for(sim in 1:no_sims){
  cat(sim, "-----\n")
  curr_sim_result <- tibble()
  
  n <- 100
  Time <- 100
  prop.groups <- c(0.2, 0.3)
  
  # intens1 <- c(1, 3, 8)/4
  # intens2 <- c(2, 3, 6)/6
  # intens <- matrix(c(intens1,intens2,intens1,intens2), 4, 3)
  intens <- matrix(c(0.25, 0.5, 1, 0.75, 1, 1, 0.50, 0.25, 1,
                     1, 0.75, 1), nrow = 4, byrow = TRUE)
  
  # intens[2, 3] <- 0.5
  # intens[3, 2] <- 1

  dynppsbm <- generateDynppsbmConst(intens,
                                    Time,
                                    n,
                                    prop.groups,
                                    directed = TRUE)
  
  ###
  # hist(dynppsbm$data$time.seq)
  proc_sim <- format_sims(sim_data = dynppsbm, n = n,
                          directed = TRUE)
  
  
  
  ### Fit SBM to count matrix
  ### this doesn't seem to be working correctly for undirected
  
  A_test <- proc_sim$edge
  events <- proc_sim$events
  
  ### Fit SBM to count matrix ###
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
  z_true <- apply(dynppsbm$z, 2, which.max)
  (sc_ari <- aricode::ARI(z_true, sc_est))
  ### this works
  
  possible_windows <- seq(from = 1, to = 10, by = 1)
  for(wind in possible_windows){
    print(wind)
    print("-----")
    ### fit Pensky Zhang to this data
    A_mats <- event_to_mat_seq(proc_sim$events,
                               Total_time = Time,
                               window_size = wind,
                               n = n)
    ## convert this to binary
    A_mats[A_mats > 0] <- 1
    
    l <- 4
    m_pz <- 4
    m0 <- 1
    l0 <- 2
    
    A_pz <- pz_estimator_3(A = A_mats, time = Time,
                           l0 = l0, m0 = m0,
                           m = m_pz)
    ### need to modify how this deals with changing the window
    ### size here, maybe don't specify r at all?
    
    # A_pz
    est_labels <- specc(A_pz, centers = 2)
    (pz_ari <- aricode::ARI(z_true, est_labels@.Data))
    
    #### then our methods, for a given window size
    #### Fit our method to this data
    K <- 2
    m <- n
    B <- matrix(runif(K * K), K, K)
    dT <- wind
    inter_T <- 1
    # capture output to not print out
    results_online <- estimate_Poisson(full_data = events,
                                       A = A_test,
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
    results_online_hawkes <- online_estimator_eff_revised(alltimes = events, 
                                                          A = A_test,
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
                                                               A = A_test,
                                                               m,
                                                               K,
                                                               H,
                                                               window,
                                                               Time,
                                                               dT,
                                                               gravity = 0.001,
                                                               MuA, tau))
    
    
    (ari_in_pois <- aricode::ARI(apply(results_online_inpois$tau, 1, which.max),
                                 z_true))
    
    ## Fitting inhomogeneous Hawkes
    
    window <- 1/2
    K <- 2 # 4 for email, 2 for college, 3 for math
    H <- 2
    #### this was the bug#####
    dT <- wind # 2 for email, 0.5 for college, 6 for math
    MuA_start <- array(runif(K * K * H), c(K, K, H))
    tau_start <- matrix(1/K, m, K)
    B_start <- matrix(runif(K * K), nrow = K, ncol = K)
    
    result_inHaw <- nonhomoHak_estimator_eff_revised(alltimes = events,
                                                     A = A_test,
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
    
    
    (ari_in_haw <- aricode::ARI(z_true, apply(result_inHaw$tau, 1, which.max)))
    
    ### then store all this
    curr_sim_result <- tibble(ARI = c(pz_ari,
                                      clust_est,
                                      haw_clust_est,
                                      ari_in_pois,
                                      ari_in_haw),
                              Method = c("PZ", "Poisson",
                                         "Hawkes", "InPois", "InHawkes"),
                              window_size = wind)
    
    all_results <- all_results %>% bind_rows(curr_sim_result)
  }
  
  ### add in the count ARI just here
  count_sim <- tibble(ARI = sc_ari,
                      Method = "Count",
                      window_size = NA)
  all_results <- all_results %>% bind_rows(count_sim)
}


saveRDS(all_results, file = here("Experiments",
                             "exp_results",
                             "fig_1_exp_1.RDS"))

