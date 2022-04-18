### Tidy Simulations for Figure 1, 
### illustrating the importance of temporal 
### structure in a community detection example
### April 18th 2022

### Here we just want to run InHomogeneous Poisson 
### For a single window setting, 
### using some initialization also (probably a small time period)
### To illustrate that community recovery can be done when
### we incorporate the changing intensity over time

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

library(kernlab)
library(ppsbm)
source(here("Experiments/", "utils.R"))
source(here("functions/", "df_to_adj.R"))
source(here("functions", "pensky_fcns.R"))
source(here("functions/init_fcn.R"))

print(here("Experiments", "thesis_output", "fig_1_exp_1.RDS"))

no_sims <- 50
Time <- 100
all_results <- tibble()
K <- 2
dT <- 1


for(sim in 1:no_sims){
  cat(sim, "-----\n")
  curr_sim_result <- tibble()
  
  n <- 100
  m <- n
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
  
  sc <- specc(A_mat, centers = K)
  sc_est <- sc@.Data
  z_true <- apply(dynppsbm$z, 2, which.max)
  (sc_ari <- aricode::ARI(z_true, sc_est))
  ### this works
  
  possible_windows <- seq(from = 1, to = 10, by = 1)
  for(wind in possible_windows){
    # print(wind)
    # print("-----")
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
    ## this seems to struggle some times
    ## put a try here?
    a <- "Error"
    attempt <- 0
    while(a == "Error" & attempt < 10){
      ## the warning here is when specc actually works
      ## so can ignore it
      a <- tryCatch(error = function(cnd) "Error",
                    specc(A_pz, centers = K)
      )
      attempt <- attempt + 1
    }
    if(attempt == 10){
      pz_ari <- NA
    } else{
      est_labels <- a ##specc(A_pz, centers = K)
      ## so don't have to run it again
      (pz_ari <- aricode::ARI(z_true, est_labels@.Data))
    }
    
    ### then store all this
    curr_sim_result <- tibble(ARI = pz_ari,
                              Method = "PZ",
                              window_size = wind,
                              sim = sim)
    
    all_results <- all_results %>% bind_rows(curr_sim_result)
  }
  ### Add the init Procedure here (using maybe just the community assignments)
  curr_n0 <- 5
  result <- dense_poisson(alltimes = events, K, n0 = curr_n0, m)
  
  ## Fitting Inhomogeneous Poisson
  dT <- 1
  H <- 2
  MuA <- array(runif(K * K * H), c(K, K, H))
  window <- 0.5
  ### construct tau
  init_tau <- matrix(0, nrow = m, ncol = K)
  for(i in seq_along(result$est_clust)){
    init_tau[i, result$est_clust[i]] <- 1
  }
  invisible(results_online_inpois <- nonhomoPois_est_init(alltimes = events,
                                                          A = A_test,
                                                          m,
                                                          K,
                                                          H,
                                                          window,
                                                          Time,
                                                          dT,
                                                          gravity = 0.001,
                                                          MuA, 
                                                          init_tau,
                                                          start = result$cut_off))
  
  
  (ari_in_pois <- aricode::ARI(apply(results_online_inpois$tau, 1, which.max),
                               z_true))
  
  ### add in the count ARI just here
  count_sim <- tibble(ARI = c(sc_ari, ari_in_pois),
                      Method = c("Count", "Point Process"),
                      window_size = NA,
                      sim = sim)
  
  all_results <- all_results %>% bind_rows(count_sim)
}

sim_settings <- list(no_sims = no_sims, 
                     Time = Time,
                     K = K,
                     prop.groups = prop.groups,
                     intens = intens,
                     windows = possible_windows)

saveRDS(sim_settings, 
        file = here("Experiments",
                    "thesis_output",
                    "setting_1_exp_1_apr_18.RDS"))

saveRDS(all_results, file = here("Experiments",
                                 "thesis_output",
                                 "fig_1_exp_1_apr_18.RDS"))

