##### Initialization Function ####
options(dplyr.summarise.inform = FALSE)
### hide the message from summarise

dense_poisson <- function(alltimes, K, n0, m) {
  
  ### select n0
  ### find node with largest degree, and its neighbours
  ###
  
  tidy_events <- as_tibble(alltimes) %>% 
    rename(send = V1, rec = V2, time = V3)
  N <- nrow(tidy_events)
  C <- 1
  # n0 <- C * log(N)
  
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
    mutate(est_lam = n/n0) %>% 
    select(rec, est_lam)
  ## check here if len(est_int) < K
  
  est <- kmeans(est_int$est_lam, centers = K)
  clust_ests <- tibble(node = est_int$rec, clust = est$cluster)
  
  neighbours <- init_events %>% 
    filter(send == max_degree) %>% 
    ### check this giving a warning too
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
      curr_est <- neighbour_events %>% 
        filter(send_clust == k1) %>% 
        filter(rec_clust == k2) %>% 
        ungroup() %>% 
        slice_sample(n = 1) %>% 
        mutate(rate = n/n0) %>% 
        pull(rate)
      ### check if empty here
      if(identical(curr_est, numeric(0))){
        init_B[k1, k2] <- 0
      }
      else{
        init_B[k1, k2] <- curr_est 
      }
    }
  }
  print(init_B)
  sorted_B <- t(apply(init_B, 1, sort))
  init_group <- rep(NA, m)
  
  for(i in (1:m)-1){
    ## 
    curr_neigh <- init_events %>% 
      filter(send == i) %>% 
      group_by(send, rec) %>% 
      summarise(count = n()) %>% 
      mutate(rate = count/n0) %>% 
      pull(rate)
    
    # then k means on these estimates
    ## check number of neighbours here too
    curr_est <- kmeans(curr_neigh, centers = K)
    curr_center <- sort(as.vector(curr_est$centers))
    
    ## then see which row of init_B this is closest to
    dists <- apply(sorted_B, 1, function(x) 
      dist(rbind(x, curr_center)))
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
        summarise(est_rate = sum(num_events)/ (n()*n0)) %>% 
        pull(est_rate)
      
      if(identical(new_est, numeric(0))){
        updated_B[k1, k2] <- 0
      }
      else{
        updated_B[k1, k2] <- new_est 
      }
      
      # updated_B[k1, k2] <- new_est
      # ### debug this
      # if(is.na(new_est)) {
      #   print("You should check the data")
      #   return(alltimes)
      # }
    }
  }
  # cat(updated_B, "\n------\n")
  print(updated_B)
  return(list(est_clust = init_group,
              est_B = updated_B,
              rest_events = remaining_events,
              cut_off = n0))
}


sparse_poisson <- function(alltimes, K, n0, m, m0){
  ###
  ###
  ### specify m0
  ###
  ### first get the event data to be used
  ###
  
  tidy_events <- as_tibble(alltimes) %>% 
    rename(send = V1, rec = V2, time = V3)
  N <- nrow(tidy_events)
  C <- 1
  # n0 <- C * log(N)
  init_events <- tidy_events %>% 
    filter(time <= n0)
  remaining_events <- alltimes[alltimes[,3]> n0, ]
  
  out_events <- init_events %>%
    group_by(send) %>%
    count() %>%
    ungroup() %>%
    rename(out = n)
  
  top_nodes <- out_events %>%
    slice_max(out, n = m0, with_ties = FALSE) %>%
    pull(send)
  
  center_matrix <- matrix(NA, nrow = m0, ncol = K)
  ### iterate along m0 nodes
  # print("Getting this far")
  for(i in seq_along(top_nodes)){
    neighs <- init_events %>%
      filter(send == top_nodes[i]) %>%
      group_by(rec) %>%
      count() %>%
      mutate(est_rate = n/n0)
    lam_neigh <- neighs %>%
      pull(est_rate)
    curr_est <- kmeans(lam_neigh, centers = K)
    ### will there be permutation issues here also?
    center_matrix[i, ] <- sort(as.vector(curr_est$centers))
  }
  # print("Here")
  ### then k means on center matrix
  init_clust <- kmeans(center_matrix, centers = K)$cluster
  
  first_clust <- tibble(node = top_nodes,
                        clust = init_clust)
  
  ### then get the events these nodes are involved in
  ### should do this the same as previous method
  init_clust_events <- init_events %>% 
    filter(send %in% top_nodes) %>% 
    filter(rec %in% top_nodes) %>% 
    group_by(send, rec) %>% 
    count() %>% 
    left_join(first_clust, by = c("send" = "node")) %>% 
    rename(send_clust = clust) %>% 
    left_join(first_clust, by = c("rec" = "node")) %>% 
    rename(rec_clust = clust)
  
  init_B <- matrix(NA, nrow = K, ncol = K)
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      curr_est <- init_clust_events %>% 
        filter(send_clust == k1) %>% 
        filter(rec_clust == k2) %>% 
        ungroup() %>% 
        slice_sample(n = 1) %>% 
        mutate(rate = n/n0) %>% 
        pull(rate)
      
      init_B[k1, k2] <- curr_est
    }
  }
  ### take all nodes not in top_nodes
  rem_nodes <- tibble(
    node = 0:(m-1),
  ) %>% 
    left_join(first_clust, by = "node") %>% 
    filter(is.na(clust)) %>% 
    pull(node)
  
  sorted_B <- t(apply(init_B, 1, sort))
  rem_group <- rep(NA, length(rem_nodes))
  
  for(i in seq_along(rem_nodes)) {
    curr_neigh <- init_events %>% 
      filter(send == rem_nodes[i]) %>% 
      group_by(send, rec) %>% 
      summarise(count = n()) %>% 
      mutate(rate = count/n0) %>% 
      pull(rate)
    
    if(length(unique(curr_neigh)) < K){
      ### randomly assign to a cluster, if number distinct values < K
      rem_group[i] <- sample(1:K, size = 1)
    }
    else if(length(unique(curr_neigh)) == K){
      ## use a different algorithm in this case
      curr_est <- kmeans(curr_neigh, centers = K, algorithm = "Lloyd")
      curr_center <- sort(as.vector(curr_est$centers))
      dists <- apply(sorted_B, 1, function(x) 
        dist(rbind(x, curr_center)))
      rem_group[i] <- which.min(dists)
    }
    else{
      curr_est <- kmeans(curr_neigh, centers = K)
      curr_center <- sort(as.vector(curr_est$centers))
      dists <- apply(sorted_B, 1, function(x) 
        dist(rbind(x, curr_center)))
      rem_group[i] <- which.min(dists)
    }
  }
  ### then update estimate of B
  ### using full clusters
  rem_clust <- tibble(node = rem_nodes,
                      clust = rem_group)
  all_clust <- bind_rows(first_clust, rem_clust)
  
  init_group <- all_clust %>% 
    arrange(node) %>% 
    pull(clust)
  
  ### the update the estimate of B
  clust_events <- init_events %>% 
    left_join(all_clust, by = c("send" = "node")) %>% 
    rename(send_clust = clust) %>% 
    left_join(all_clust, by = c("rec" = "node")) %>% 
    rename(rec_clust = clust) 
  
  updated_B <- init_B
  for(k1 in 1:K){
    for(k2 in 1:K){
      new_est <- clust_events %>% 
        filter(send_clust == k1) %>% 
        filter(rec_clust == k2) %>% 
        group_by(send, rec) %>% 
        summarise(num_events = n()) %>% 
        ungroup() %>% 
        ### fitting a common Poisson to each of these pairs
        summarise(est_rate = sum(num_events)/ (n()*n0)) %>% 
        pull(est_rate)
      updated_B[k1, k2] <- new_est
    }
  }
  # print(updated_B)
  
  ### then return all this stuff in a list
  return(list(est_clust = init_group,
              est_B = updated_B,
              rest_events = remaining_events,
              cut_off = n0))
}


