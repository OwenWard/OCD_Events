#### helper functions for cleaning data and loading 
#### required functions, packages

library(tidyverse)
library(Rcpp)
library(ppsbm)
library(kernlab)
library(RcppArmadillo)
sourceCpp(here("cpp_files", "minimal_cpp","minimal_functions.cpp"))
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
  tidy_sim <- as_tibble(sim_trip, .name_repair = ~ c("V1", "V2", "V3")) %>% 
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


#' Generate Network Point Process Data
#' 
#' Modified version of function `ppsbm::generateDynppsbmConst`, 
#' allowing sparsity in the edges which interact
#'
#' @param intens 
#' @param Time 
#' @param n 
#' @param prop.groups 
#' @param directed 
#' @param prob_edge 
#'
#' @return
#' @export
#'
#' @examples
gen_ppsbm <- function(intens, Time, n, prop.groups, directed = TRUE, prob_edge = 1) 
{
  Q <- length(prop.groups)
  N_Q <- if (directed) 
    Q^2
  else Q * (Q + 1)/2
  if (nrow(intens) != N_Q) 
    stop("not a correct number of intensities")
  z <- rmultinom(n, 1, prop.groups)
  group <- Rfast::colMaxs(z, FALSE)
  time.seq <- NULL
  type.seq <- NULL
  vec_i <- if (directed) 
    1:n
  else 1:(n - 1)
  for (i in vec_i) {
    vec_j <- if (directed) 
      (1:n)[-i]
    else (i + 1):n
    for (j in vec_j) {
      edge <- rbernoulli(1, p = prob_edge)
      if(edge){
        type.ql <- convertGroupPair(group[i], group[j], Q, 
                                    directed)
        proc <- generatePPConst(intens[type.ql, ], Time)
        if (length(proc) > 0) {
          time.seq <- c(time.seq, proc)
          type.seq <- c(type.seq, rep(convertNodePair(i, 
                                                      j, n, directed), length(proc)))
        }
      }
    }
  }
  ordre <- order(time.seq)
  time.seq <- time.seq[ordre]
  type.seq <- type.seq[ordre]
  data <- list(time.seq = time.seq, type.seq = type.seq, Time = Time)
  return(list(data = data, z = z, intens = intens))
}