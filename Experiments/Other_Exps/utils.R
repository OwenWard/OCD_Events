#### helper functions for cleaning data and loading 
#### required functions, packages

library(here)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
sourceCpp(here("cpp_files", "onlineblock.cpp"))
# library(ppsbm)
theme_set(theme_bw())

#### function to process output from ppsbm

format_sims <- function(sim_data, n, directed = TRUE) {
  all_pairs <- listNodePairs(n = n, directed = directed)
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
      filter(send == i - 1) %>% # to account for zero indexing here
      distinct(rec) %>% 
      arrange(rec) %>% 
      pull(rec)
    A_new[[i]] <- edge
  }
  list(events = sim_trip, edge = A_new)
}
