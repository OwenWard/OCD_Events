#### helper functions

library(here)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
sourceCpp(here("cpp_files", "onlineblock.cpp"))
library(ppsbm)

#### function to process output from ppsbm

format_sims <- function(sim_data, n) {
  all_pairs <- listNodePairs(n = n)
  event_pairs <- all_pairs[sim_data$data$type.seq,]
  sim_trip <- cbind(event_pairs, sim_data$data$time.seq)
  sim_trip[, 1] <- sim_trip[, 1] - 1
  sim_trip[, 2] <- sim_trip[, 2] - 1
  # zero index for Rcpp
  
  # construct edge list
  tidy_sim <- as_tibble(sim_trip) %>% 
    rename(send = V1, rec = V2, time = V3)
  A_new <- list()
  for(i in 1:n){
    edge <- tidy_sim %>% 
      filter(send == i) %>% 
      distinct(rec) %>% 
      arrange(rec) %>% 
      pull(rec)
    A_new[[i]] <- edge
  }
  list(events = sim_trip, edge = A_new)
}
