##### Initialization Function ####
options(dplyr.summarise.inform = FALSE)
### hide the message from summarise

dense_poisson <- function(alltimes, K) {
  
  ### select n0
  ### find node with largest degree, and its neighbours
  ###
  
  tidy_events <- as_tibble(alltimes) %>% 
    rename(send = V1, rec = V2, time = V3)
  N <- nrow(tidy_events)
  C <- 1
  n0 <- C * log(N)
  
  init_events <- tidy_events %>% 
    filter(time <= n0)
  remaining_events <- alltimes[alltimes[,3]> n0, ]
  
  out_events <- init_events %>% 
    group_by(send) %>% 
    count() %>% 
    ungroup() %>% 
    rename(out = n)
  
  max_degree <- out_events %>% 
    slice_max(out, n = 1) %>% 
    pull(send)
  
  est_int <- init_events %>% 
    filter(send == max_degree) %>% 
    group_by(rec) %>% 
    count() %>% 
    mutate(est_lam = n/Time) %>% 
    select(rec, est_lam)
  
  est <- kmeans(est_int$est_lam, centers = K)
  clust_ests <- tibble(node = est_int$rec, clust = est$cluster)
  
  neighbours <- init_events %>% 
    filter(send == max_degree) %>% 
    distinct(rec) %>% 
    arrange(rec) %>% 
    pull(rec)
  
  ### then get all events between these pairs
  
  neighbour_events <- init_events %>% 
    filter(send %in% neighbours) %>% 
    filter(rec %in% neighbours) %>% 
    group_by(send, rec) %>% 
    count() %>% 
    ## then join the cluster assignments here
    left_join(clust_ests, by = c("send" = "node")) %>% 
    rename(send_clust = clust) %>% 
    left_join(clust_ests, by = c("rec" = "node")) %>% 
    rename(rec_clust = clust)
  
  init_B <- matrix(NA, nrow = K, ncol = K)
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      init_B[k1, k2] <- neighbour_events %>% 
        filter(send_clust == k1) %>% 
        filter(rec_clust == k2) %>% 
        ungroup() %>% 
        slice_sample(n = 1) %>% 
        mutate(rate = n/Time) %>% 
        pull(rate)
    }
  }

  init_group <- rep(NA, m)
  
  for(i in (1:m)-1){
    ## 
    curr_neigh <- init_events %>% 
      filter(send == i) %>% 
      group_by(send, rec) %>% 
      summarise(count = n()) %>% 
      mutate(rate = count/Time) %>% 
      pull(rate)
    
    # then k means on these estimates
    curr_est <- kmeans(curr_neigh, centers = K)
    curr_center <- curr_est$centers
    
    ## then see which row of init_B this is closest to
    dists <- apply(init_B, 1, function(x) 
      dist(rbind(x, as.vector(curr_center))))
    init_group[i + 1] <- which.min(dists)
    ## to deal with the zero indexing of the raw data
  }
  
  ### then the final part, the final estimate of B
  ## add clust assignment to tidy events
  
  clust_info <- tibble(nodes = 0:(m-1), clust = init_group)
  
  clust_events <- init_events %>% 
    left_join(clust_info, by = c("send" = "nodes")) %>% 
    rename(clust_send = clust) %>% 
    left_join(clust_info, by = c("rec" = "nodes")) %>% 
    rename(clust_rec = clust)
  
  updated_B <- init_B
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      new_est <- clust_events %>% 
        filter(clust_send == k1) %>% 
        filter(clust_rec == k2) %>% 
        group_by(send, rec) %>% 
        summarise(num_events = n()) %>% 
        ungroup() %>% 
        ### fitting a common Poisson to each of these pairs
        summarise(est_rate = sum(num_events)/ (n()*Time)) %>% 
        pull(est_rate)
      updated_B[k1, k2] <- new_est
    }
  }
  return(list(est_clust = init_group,
              est_B = updated_B,
              rest_events = remaining_events,
              cut_off = n0))
}


### run through and apply

K <- 2
m <- 100
Time <- 100
sparsity <- 0.5
true_Mu <- matrix(c(0.5, 0.05, 0.05, 1), 
                  nrow = 2, ncol = 2, byrow = T)
## excitation, if used (for Hawkes)
true_B <- matrix(c(0.5, 0, 0, .5), nrow = K, ncol = K, byrow = TRUE)
if(model == "Poisson") {
  true_B <- matrix(0, K, K)
}
Pi <- matrix(c(0.5, 0.5), 1, 2)
Z <- c(rep(0, m * Pi[1]), rep(1, m * Pi[2]))
# then generate A
A <- list()
for(i in 1:m){
  # could sample these here with SBM structure...
  num_edge = m * sparsity
  edge <- sample(m, num_edge) - 1
  edge <- sort(edge)
  A[[i]] <- edge
}


alltimes <- sampleBlockHak(Time, A, Z, Mu = true_Mu, B = true_B, lam = 1)


result <- dense_poisson(alltimes, K = 2)

aricode::ARI(result$est_clust, Z)

