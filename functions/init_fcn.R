##### Initialization Function ####
options(dplyr.summarise.inform = FALSE)
### hide the message from summarise



## for a dense point process network, where data believed to 
## follow block homogeneous Poisson process, provide initial 
## estimates of the community structure, rate matrix

dense_poisson <- function(alltimes, K, n0, m) {
  
  ### select n0
  ### find node with largest degree, and its neighbours
  ###
  
  tidy_events <- as_tibble(alltimes, 
                           .name_repair = ~ c("V1", "V2", "V3")) %>% 
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
  
  ## check here if len(est_int$est_lam) < K
  lam_vec <- est_int$est_lam
  # print(lam_vec)
  if(length(unique(lam_vec)) < K){
    ### give random estimates
    est <- sample(1:K, size = length(lam_vec), replace = TRUE)
  }else if(length(unique(lam_vec)) == K){
    est <- kmeans(lam_vec, centers = K, algorithm = "Lloyd")
  }else{
    est <- kmeans(lam_vec, centers = K) 
  }
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
  # print(init_B)
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
    
    ## deal with length(curr_neigh) <= K
    if(length(unique(curr_neigh)) < K){
      # print("Need to consider this case")
      curr_est <- c(rep(0, K-length(unique(curr_neigh))), 
                    unique(curr_neigh))
      curr_center <- sort(curr_est)
      ### pad with zeros to have length K and then sort?
    } else if(length(unique(curr_neigh)) == K){
      curr_est <- kmeans(curr_neigh, centers = K, algorithm = "Lloyd")
      curr_center <- sort(as.vector(curr_est$centers))
    } else{
      curr_est <- kmeans(curr_neigh, centers = K)
      curr_center <- sort(as.vector(curr_est$centers))
    }
    
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
      else if(is.nan(new_est)){
        updated_B[k1, k2] <- runif(1)
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
  # print(updated_B)
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
  
  tidy_events <- as_tibble(alltimes, 
                           .name_repair = ~ c("V1", "V2", "V3")) %>% 
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
  
  center_matrix <- matrix(NA, nrow = length(top_nodes), ncol = K)
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
    ### dealing with cases where length(lam_neigh) < K
    if(length(unique(lam_neigh)) < K){
      # print("Need to consider this case")
      curr_est <- c(rep(0, K-length(unique(lam_neigh))), 
                    unique(lam_neigh))
      center_matrix[i, ] <- sort(curr_est)
      ### pad with zeros to have length K and then sort?
    } else if(length(unique(lam_neigh)) == K){
      curr_est <- kmeans(lam_neigh, centers = K, algorithm = "Lloyd")
      center_matrix[i, ] <- sort(as.vector(curr_est$centers))
    } else{
      curr_est <- kmeans(lam_neigh, centers = K)
      center_matrix[i, ] <- sort(as.vector(curr_est$centers))
    }
    
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
      curr_est <- tryCatch(
        error = function(cnd) runif(n = 1),
        init_clust_events %>% 
          filter(send_clust == k1) %>% 
          filter(rec_clust == k2) %>% 
          ungroup() %>% 
          slice_sample(n = 1) %>% 
          mutate(rate = n/n0) %>% 
          pull(rate)
      )
      
      if(length(curr_est) > 0) {
        init_B[k1, k2] <- curr_est
      }
      else{
        init_B[k1, k2] <- 0
      }
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
  # print(updated_B)
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
      if(is.nan(new_est)){
        updated_B[k1, k2] <- runif(1)
        ## random init if doesn't work
      }else{
        updated_B[k1, k2] <- new_est 
      }
    }
  }
  # print(updated_B)
  
  ### then return all this stuff in a list
  return(list(est_clust = init_group,
              est_B = updated_B,
              rest_events = remaining_events,
              cut_off = n0))
}

## helper functions for inhom Poisson

#' Count the number of events in a given windowed point process
#'
#' @param events 
#' @param jumps 
#' @param which_h 
#' @param window 
#'
#' @return list, containing the number of events in each window, number for a
#' given basis function and the time each basis function is observed for
#' @export
#'
#' @examples
num_in_wind <- function(events, jumps, which_h, window){
  ## given the breakpoints (here right is upper limit,)
  ## return vector of length length(jumps)-1,
  ## giving how many events in a specific window 
  counts <- rep(0, length(jumps) - 1)
  time_H <- rep(0,H)
  counts_H <- rep(0, H)
  for(i in seq_along(jumps[-length(jumps)])){
    lower <- jumps[i]
    upper <- jumps[i+1]
    counts[i] <- length(which(events < upper & events > lower))
    counts_H[which_h[i]] <- counts_H[which_h[i]] + counts[i]  
    time_H[which_h[i]] <- time_H[which_h[i]] + window 
  }
  list(counts = counts, time = time_H, counts_H = counts_H)
  ## only change to make this more general
  ## would to instead have each entry be a vector of 
  ## the event times in a window rather than
  ## just the count
}

#' Estimate the rate functions of an inhomogeneous Poisson Process
#'
#' @param t_start 
#' @param t_end 
#' @param events 
#' @param window 
#' @param H 
#' 
#' @description This rounds down to windows which are fully observed 
#' based on the inputs
#'
#' @return the estimated rates
#' @export
#'
#' @examples
fit_inpois <- function(events, t_start, t_end, window, H){
  
  ## bin the data based on windows, then just concatenate
  ## it and fit separately?
  ## this won't work for hawkes because history is common,
  ## or would that just be an additional step then?
  ## return the estimated rates for each of the corresponding windows
  
  num_events <- rep(0, H)
  # obs_time <- rep(0, H)
  h1 <- floor(t_start/window) # start closest and to the left 
  h2 <- floor(t_end/window) # this will throw out incomplete intervals
  jumps <- h1:h2 * window # to get the actual jumps
  which_h <- rep(1:H, length.out = length(jumps) - 1)
  ## then bin the data based on these jumps
  count_data <- num_in_wind(events, jumps, which_h, window)
  obs_time <- count_data$time
  H_counts <- count_data$counts_H
  est_rates <- H_counts / obs_time
  est_rates
}


dense_inhom_Poisson <- function(alltimes, K, H, window, t_start, n0, m) {
  
  ## given a choice of H, window,
  ## estimate the individual rates for each node in each window
  ## then do K-mean on those H-dimensional vectors
  ## to get cluster assignments, corresponding updated estimates 
  ## of rates
  
  tidy_events <- as_tibble(alltimes,
                           .name_repair = ~ c("V1", "V2", "V3")) %>% 
    rename(send = V1, rec = V2, time = V3)
  N <- nrow(tidy_events)
  C <- 1
  
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
  
  ## need to change this below
  ## want to fit based on the window
  
  max_events <- init_events %>% 
    filter(send == max_degree) %>% 
    group_by(rec) %>%
    summarise(events = list(time)) 
  names(max_events$events) <- max_events$rec
  ests <- map_dfr(max_events$events, fit_inpois, t_start = t_start,
                  t_end = n0, window, H)
  est_names <- paste0("H", 1:H)
  ests <- as_tibble(t(ests), .name_repair = ~ c(paste0("V", 1:H)))
  ## then do clustering here
  
  ## below pulled from hom poisson dense
  if(nrow(unique(ests)) < K){
    ### give random estimates
    est <- sample(1:K, size = nrow(ests), replace = TRUE)
  }else if(length(unique(ests)) == K){
    est <- kmeans(ests, centers = K, algorithm = "Lloyd")
  }else{
    est <- kmeans(ests, centers = K) 
  }
  clust_ests <- tibble(node = max_events$rec, clust = est$cluster)
  
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
    summarise(events = list(time)) %>% 
    ## then join the cluster assignments here
    left_join(clust_ests, by = c("send" = "node")) %>% 
    rename(send_clust = clust) %>% 
    left_join(clust_ests, by = c("rec" = "node")) %>% 
    rename(rec_clust = clust)
  
  init_MuA <- array(NA, dim = c(K, K, H))
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      curr_data <- neighbour_events %>% 
        filter(send_clust == k1) %>% 
        filter(rec_clust == k2) %>% 
        ungroup() 
      if(nrow(curr_data) > 0){
        curr_est <- curr_data %>% 
          slice_sample(n = 1) %>%
          ## change this bit here based on the function being used
          mutate(ests = list(fit_inpois(events[[1]], t_start = t_start,
                                        t_end = n0, window, H))) %>% 
          pull(ests)
          
      }
      else{
        curr_est <- runif(H)
      }
        ## deal with empty here also
      ### check if empty here
      if(identical(unlist(curr_est), numeric(0))){
        ## update this to update the array entry
        init_MuA[k1, k2, ] <- 0
      }
      else{
        ## same here
        init_MuA[k1, k2, ] <- unlist(curr_est) 
      }
    }
  }
  ord <- order(init_MuA[1, , 1])
  sorted_MuA <- apply(init_MuA, c(1, 3), function(x) x[ord])
  # works but bad example here because same values 
  ## could find first which not tied?
  init_group <- rep(NA, m)
  
  ## stack into K *(H*K) matrix
  Mu_matrix <- matrix(as.vector(aperm(sorted_MuA, c(3, 2, 1))),
                      nrow = 2, ncol = 4, byrow = T)
  
  
  for(i in (1:m)-1){
    curr_neigh <- init_events %>% 
      filter(send == i) %>% 
      group_by(send, rec) %>% 
      summarise(events = list(time)) 
    ## same as above
    names(curr_neigh$events) <- curr_neigh$rec
    ests <- map_dfr(curr_neigh$events, fit_inpois, t_start = t_start,
                    t_end = n0, window, H)
    ests <- as_tibble(t(ests), .name_repair = ~ c(paste0("V", 1:H)))
    
    ## TO DO do I need to sort in here?
    ## deal with length(curr_neigh) <= K
    if(nrow(unique(ests)) < K){
      # print("Need to consider this case")
      curr_est <- c(rep(0, K-length(unique(ests))), 
                    unique(ests))
      curr_center <- as.vector(t(curr_est))
      ### pad with zeros to have length K and then sort?
    } else if(nrow(unique(ests)) == K){
      curr_est <- kmeans(ests, centers = K, algorithm = "Lloyd")
      curr_center <- as.vector(t(curr_est$centers))
    } else{
      curr_est <- kmeans(ests, centers = K)
      ## here each row is a basis function, want to 
      ## stretch this out to a vector, putting rows side by side
      curr_center <- as.vector(t(curr_est$centers))
    }
    
    ## then see which row of init_B this is closest to
    dists <- apply(Mu_matrix, 1, function(x) 
      dist(rbind(x, curr_center)))
    init_group[i + 1] <- which.min(dists)
    ## to deal with the zero indexing of the raw data
    ## dist from init_MuA[i, j, ] to curr_center
  }
  
  ### then the final part, the final estimate of B
  ## add clust assignment to tidy events
  
  clust_info <- tibble(nodes = 0:(m-1), clust = init_group)
  
  clust_events <- init_events %>% 
    left_join(clust_info, by = c("send" = "nodes")) %>% 
    rename(clust_send = clust) %>% 
    left_join(clust_info, by = c("rec" = "nodes")) %>% 
    rename(clust_rec = clust)
  
  updated_Mu <- sorted_MuA
  t_end <- n0
  h1 <- floor(t_start/window) # start closest and to the left 
  h2 <- floor(t_end/window) # this will throw out incomplete intervals
  jumps <- h1:h2 * window # to get the actual jumps
  which_h <- rep(1:H, length.out = length(jumps) - 1)
  
  ## occasional bug here, think if one group empty?
  print(table(init_group))
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      curr_data <- clust_events %>% 
        filter(clust_send == k1) %>% 
        filter(clust_rec == k2) %>% 
        group_by(send, rec) %>% 
        summarise(events = list(time)) %>%
        rowwise() %>% 
        mutate(counts = list(num_in_wind(unlist(events),
                                         jumps, which_h, window))) %>% 
        ungroup() 
      
      if(nrow(curr_data) > 0) {
        curr_data %>% 
          unnest_wider(col = counts) %>% 
          unnest_wider(col = c("counts_H", "time"), names_sep = "_")
        ## then sum across counts and divide by all the time
        total_counts <- curr_data %>% select(starts_with("counts_H")) %>% 
          summarise_all(sum) %>% 
          as.numeric()
        
        total_time <- curr_data %>% select(starts_with("time_")) %>% 
          summarise_all(sum) %>% 
          as.numeric()
        new_est <- total_counts/total_time
      }
      else{
        new_est <- runif(H)
      }
      
      if(identical(new_est, numeric(0))){
        updated_Mu[k1, k2, ] <- 0
      }
      # else if(is.nan(new_est)){
      #   updated_Mu[k1, k2, ] <- runif(1)
      # }
      else{
        updated_Mu[k1, k2, ] <- new_est 
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
  # print(updated_B)
  return(list(est_clust = init_group,
              est_Mu = updated_Mu,
              rest_events = remaining_events,
              cut_off = n0))
  
}
