negloglik_hp_comm <- function(vec, events, end = end) {
  negllh <- 0
  # transforms input list object into vector so that it can be used in optim
  for(i in seq_along(events)){
    curr_events <- unlist(events[i])
    object <- list(lambda0 = vec[1], alpha = vec[2], beta = vec[3])
    class(object) <- "hp"
    negllh <- negllh + ppdiag:::negloglik.hp(object = object,
                                             events = curr_events, end = end)
  }
  return(negllh)
}


fithp_common <- function(events, end = end, vec = c(0.1, 0.2, 0.3)) {
  con_mat <- matrix(0, nrow = 4, ncol = 3)
  con_mat[1:3, 1:3] <- diag(1, nrow = 3) # constraint all positive
  con_mat[4, ] <- c(0, -1, 1) # constrain beta > alpha
  const <- rep(0, 4)
  start <- 0
  if (start == end) {
    stop("Start and end time are equal.")
  }
  hawkes.par <- constrOptim(vec,
                            f = negloglik_hp_comm,
                            events = events, end = end,
                            ui = con_mat,
                            ci = const,
                            method = "Nelder-Mead"
  )
  a <- hawkes.par$par[2]
  b <- hawkes.par$par[3]
  c <- hawkes.par$convergence[1]
  i <- 1
  while (a >= b) {
    if (c == 1 | c == 10) {
      # if maxit reached or degeneracy, then we automatically refit
      hawkes.par <- constrOptim(
        par = vec, fn = negloglik_hp_comm,
        events = events, end = max(events),
        control = list(maxit = 1000),
        ui = con_mat,
        ci = const,
        method = "Nelder-Mead"
      )
      a <- hawkes.par$par[2]
      b <- hawkes.par$par[3]
      c <- hawkes.par$convergence[1]
    }
    else {
      # if maxit not reached but a>=b,
      # then we refit
      hawkes.par <- constrOptim(
        par = vec, fn = negloglik_hp_comm,
        events = events, end = max(events),
        control = list(maxit = 1000),
        ui = con_mat,
        ci = const,
        method = "Nelder-Mead"
      )
      a <- hawkes.par$par[2]
      b <- hawkes.par$par[3]
      c <- hawkes.par$convergence[1]
    }
    i <- i + 1
    if (i > 10) {
      stop("No solution after 10 attempts. Try a different initial vector")
    }
  }
  return(list(lambda0 = hawkes.par$par[1],
              alpha = hawkes.par$par[2],
              beta = hawkes.par$par[3]))
}






fithp_init <- function(events, end = max(events), vec = c(0.1, 0.2, 0.3)) {
  con_mat <- matrix(0, nrow = 4, ncol = 3)
  con_mat[1:3, 1:3] <- diag(1, nrow = 3) # constraint all positive
  con_mat[4, ] <- c(0, -1, 1) # constrain beta > alpha
  const <- rep(0, 4)
  start <- 0
  if (start == end) {
    stop("Start and end time are equal.")
  }
  hawkes.par <- constrOptim(vec,
                            f = ppdiag:::negloglik_hp,
                            events = events, end = end,
                            ui = con_mat,
                            ci = const,
                            method = "Nelder-Mead"
  )
  a <- hawkes.par$par[2]
  b <- hawkes.par$par[3]
  c <- hawkes.par$convergence[1]
  i <- 1
  while (a >= b) {
    if (c == 1 | c == 10) {
      # if maxit reached or degeneracy, then we automatically refit
      hawkes.par <- constrOptim(
        par = vec, fn = ppdiag:::negloglik_hp,
        events = events, end = max(events),
        control = list(maxit = 1000),
        ui = con_mat,
        ci = const,
        method = "Nelder-Mead"
      )
      a <- hawkes.par$par[2]
      b <- hawkes.par$par[3]
      c <- hawkes.par$convergence[1]
    }
    else {
      # if maxit not reached but a>=b,
      # then we refit
      hawkes.par <- constrOptim(
        par = vec, fn = ppdiag:::negloglik_hp,
        events = events, end = max(events),
        control = list(maxit = 1000),
        ui = con_mat,
        ci = const,
        method = "Nelder-Mead"
      )
      a <- hawkes.par$par[2]
      b <- hawkes.par$par[3]
      c <- hawkes.par$convergence[1]
    }
    i <- i + 1
    if (i > 10) {
      stop("No solution after 10 attempts. Try a different initial vector")
    }
  }
  return(list(lambda0 = hawkes.par$par[1],
              alpha = hawkes.par$par[2],
              beta = hawkes.par$par[3]))
}



sparse_Hawkes <- function(alltimes, K, n0, m, m0){
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
  
  center_matrix <- matrix(NA, nrow = length(top_nodes), ncol = K * 3)
  ### iterate along m0 nodes
  # print("Getting this far")
  for(i in seq_along(top_nodes)){
    ## need to modify here to fit Hawkes rather than Poisson
    neighs <- init_events %>%
      filter(send == top_nodes[i]) %>%
      group_by(rec) %>%
      summarise(events = list(time)) %>% 
      rowwise() %>% 
      ## CHECK THIS
      mutate(est_fit = list(fithp_init(events = unlist(events), end = n0))) %>% 
      unnest_wider(col = est_fit)
      ## updating here, need to figure out how to get the 
      ## ests correctly, should be very similar to inhom Poisson though
    lam_neigh <- neighs %>%
      select(lambda0:beta)
    ### dealing with cases where length(lam_neigh) < K
    if(nrow(unique(lam_neigh)) < K){
      # print("Need to consider this case")
      curr_est <- c(rep(0, 3 * K - 3 *nrow(unique(lam_neigh))), 
                    as.vector(t(unique(lam_neigh))))
      center_matrix[i, ] <- curr_est
      ### pad with zeros to have length K and then sort?
    } else if(nrow(unique(lam_neigh)) == K){
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
  
  ## update this to an array, for each of the parameters to be estimated
  ## not assuming common beta here
  init_B <- array(NA, dim = c(K, K, 3))
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      curr_est <- tryCatch(
        error = function(cnd) runif(n = 3),
        init_clust_events %>% 
          filter(send_clust == k1) %>% 
          filter(rec_clust == k2) %>% 
          ungroup() %>% 
          slice_sample(n = 1) %>% 
          ### update here to hawkes also
          mutate(est_fit = list(fithp_init(events = unlist(events), end = n0))) %>% 
          unnest_wider(col = est_fit) %>% 
          select(lambda0: beta) %>% 
          as.numeric()
      )
      
      if(length(curr_est) > 0) {
        init_B[k1, k2, ] <- curr_est
      }
      else{
        init_B[k1, k2, ] <- 0
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
  
  ## redo the ordering here, in terms of the base rates (for example)
  ord <- order(init_B[1, , 1])
  sorted_B <- apply(init_B, c(1, 3), function(x) x[ord])
  
  B_matrix <- matrix(as.vector(aperm(sorted_B,
                                     c(3, 2, 1))),
                     nrow = K, ncol = K * 3, byrow = TRUE)
  rem_group <- rep(NA, length(rem_nodes))
  
  for(i in seq_along(rem_nodes)) {
    curr_neigh <- init_events %>% 
      filter(send == rem_nodes[i]) %>% 
      group_by(send, rec) %>% 
      summarise(events = list(time)) %>%
      rowwise() %>% 
      mutate(est_fit = list(fithp_init(events = unlist(events),
                                       end = n0))) %>% 
      unnest_wider(col = est_fit) %>% 
      select(lambda0:beta)
    
    if(nrow(unique(curr_neigh)) < K){
      ### randomly assign to a cluster, if number distinct values < K
      rem_group[i] <- sample(1:K, size = 1)
    }
    else if(nrow(unique(curr_neigh)) == K){
      ## use a different algorithm in this case
      curr_est <- kmeans(curr_neigh, centers = K, algorithm = "Lloyd")
      curr_center <- sort(as.vector(curr_est$centers))
      dists <- apply(B_matrix, 1, function(x) 
        dist(rbind(x, curr_center)))
      rem_group[i] <- which.min(dists)
    }
    else{
      curr_est <- kmeans(curr_neigh, centers = K)
      curr_center <- sort(as.vector(curr_est$centers))
      dists <- apply(B_matrix, 1, function(x) 
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
        summarise(events = list(time))
        ## now need to fit a common hawkes here, which will take 
        ## a bit more work
        new_lam <- as.numeric(fithp_common(new_est$events, end = n0))
      if(identical(new_lam, c(0.1, 0.2, 0.3))){
        updated_B[k1, k2, ] <- runif(3)
        ## random init if doesn't work
      }else{
        updated_B[k1, k2, ] <- new_lam 
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
