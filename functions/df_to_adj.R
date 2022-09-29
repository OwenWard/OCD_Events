#### count tibble to Adjacency matrix

summ_to_adj_mat <- function(df, n = NA) {
  ### assume send and rec are zero indexed in df, but obviously won't be in
  send <- df$Send
  rec <- df$Rec
  nodes <- unique(sort(c(send, rec)))
  if(!is.na(n)) {
    adj <- matrix(0, nrow = n, ncol = n)
  }
  else{
    adj <- matrix(0, nrow = length(nodes), ncol = length(nodes))
    
  }
  ### then just a for loop over df
  for(i in 1:nrow(df)){
    curr_send <- df$Send[i]
    curr_rec <- df$Rec[i]
    adj[curr_send + 1, curr_rec + 1 ] <- df$n[i]
  }
  adj
}



#' Title
#'
#' @param alltimes the raw data
#' @param Total_time the total observation period
#' @param window_size how big each observation window should be
#'
#' @return A_seq, an array of adjacency matrices
#' @export
#'
#' @examples
event_to_mat_seq <- function(alltimes, Total_time, window_size, n) {
  ### assume that Total_time/window goes nicely
  
  tidy_events <- as_tibble(alltimes) %>% 
    rename(Send = V1,
           Rec = V2,
           Time = V3)
      
  n_mats <- Total_time/window_size
  A_list <- array(0, dim = c(n, n, n_mats))
  end_point <- seq(from = window_size, to = Total_time, by = window_size)
  for(i in seq_along(end_point)){
    end_time <- end_point[i]
    start_time <- end_time - window_size
    curr_events <- tidy_events %>% 
      filter(Time > start_time) %>% 
      filter(Time < end_time) %>% 
      group_by(Send, Rec) %>% 
      count()
    curr_A <- summ_to_adj_mat(curr_events, n)
    A_list[, , i] <- curr_A
    
  }
  A_list
}