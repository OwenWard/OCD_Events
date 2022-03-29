###### Functions needed for estimators in Pensky paper ####
## This script contains the code needed to construct the
## complex estimators needed in their algorithm

## First section to estimate the weights W of each
## observation

# this function appears correct
P_coefs <- function(h, k, r) {
  value <- 0
  for (i in 0:(r - 1)) {
    value <- value + (i^h) * ((r - i)^k)
  }
  value <- value * (r^(-h - k - 1))
  return(value)
}

K_matrix_1 <- function(l0, m0, m, r) {
  output <- matrix(0, nrow = l0 + 1, ncol = m + 1)
  # entries in the first row
  for (j in 1:(m + 1)) {
    output[1, j] <- 1 + 2 * r * P_coefs(h = j - 1 + m0, k = 0, r)
  }
  for (i in 2:(l0 + 1)) {
    for (j in 1:(m + 1)) {
      output[i, j] <- P_coefs(h = j - 1 + m0, k = 2 * (i - 1), r)
      # seeing as a goes a_0 to a_m
    }
  }
  return(output)
}


K_matrix_3 <- function(l0, m0, m, r) {
  output <- matrix(0, nrow = l + 1, ncol = m + 1)
  # entries in the first row
  for (j in 1:(m + 1)) {
    output[1, j] <- 1 + r * P_coefs(h = j - 1 + m0, k = 0, r)
  }
  for (i in 2:(l + 1)) {
    for (j in 1:(m + 1)) {
      output[i, j] <- P_coefs(h = j - 1 + m0, k = i - 1, r)
      # seeing as a goes a_0 to a_m
    }
  }
  return(output)
}


## function to get the weights on a sequence of
## adjacency matrices given s,r,t,T,K
## these are currently only for D_1

get_a_vec_1 <- function(l0, m0, m, r) {
  # this returns the vector a with m+1 elements
  K <- K_matrix_1(l0, m0, m, r)
  e <- rep(0, m + 1)
  e[1] <- 2 * r + 1
  result <- solve(K) %*% e
  return(result)
}

get_a_vec_3 <- function(l0, m0, m, r) {
  # this returns the vector a with m+1 elements
  K <- K_matrix_3(l0, m0, m, r)
  e <- rep(0, m + 1)
  e[1] <- r + 1
  result <- solve(K) %*% e
  return(result)
}

get_W_1 <- function(a_vec, r, m0) {
  output <- c()
  for (i in -r:r) {
    val <- 0
    for (j in 0:m) {
      temp <- a_vec[j + 1] * (r - abs(i))^(j + m0) * r^(-j - m0)
      val <- val + temp
    }
    output <- c(output, val)
  }
  return(output)
}

get_W_3 <- function(a_vec, r, m0, m) {
  output <- c()
  for (i in (-r):0) {
    val <- 0
    for (j in 0:m) {
      temp <- a_vec[j + 1] * (r + i)^(j + m0) * r^(-j - m0)
      val <- val + temp
    }
    output <- c(output, val)
  }
  return(output)
}

## now a function to perform spectral clustering
pz_estimator_1 <- function(A, time, l0, m0, m, r) {
  # takes in array of adjacency matrices and constructs
  # the estimator at a given time with window of size r.

  # first construct the coefficients
  a <- get_a_vec_1(l0, m0, m, r)
  w <- get_W_1(a, r, m0)

  # then extract the corresponding part of A
  n <- dim(A)[1]
  A_wind <- A[, , (time - r):(time + r)]
  A_agg <- matrix(0, nrow = n, ncol = n)
  num_winds <- length(w)
  for (i in 1:num_winds) {
    A_agg <- A_agg + w[i] * A_wind[, , i]
  }
  return(A_agg / num_winds)
}

pz_estimator_3 <- function(A, time, l0, m0, m) {
  # takes in array of adjacency matrices and constructs
  # the estimator at a given time with window of size r.

  r <- dim(A)[3] - 1
  ## this allows it to be more general, match result when specifying 
  ## r explicitly 
  ## always use all data
  # first construct the coefficients
  a <- get_a_vec_3(l0, m0, m, r)
  w <- get_W_3(a, r, m0, m)

  # then extract the corresponding part of A
  n <- dim(A)[1]
  A_wind <- A[, , ]
  ## enforce the use of all entries in A
  A_agg <- matrix(0, nrow = n, ncol = n)
  num_winds <- length(w)
  for (i in 1:num_winds) {
    A_agg <- A_agg + w[i] * A_wind[, , i]
  }
  return(A_agg / num_winds)
}


##
## to be completed ####
# K_matrix_3 <- function(l,m0,m,r){
#   # what equality do we need here to
#   # eunsure a unique solution?
#   output<- matrix(0,nrow = l+1,ncol = m+1)
#   for(j in 1:(m+1)){
#     # these have to be updated
#     output[1,j] <- 1+r*P_coefs(h=j+m0,k=0,r)
#   }
#   for(i in 2:(l+1)){
#     for(j in 1:(m+1)){
#       # these have to be updated
#       output[i,j] <- P_coefs(h=j+m0,k=i,r)
#     }
#   }
#   return(output)
# }

## debugging code ####
# l0 <- 1
# l <- 2*l0
# m <- l # need m=l0 for unique solution, not in general
# m0 <- 1
# r <- 6
#
# #K <- K_matrix_1(l0,m0,m,r)
# K <- K_matrix_3(l0,m0,m,r)
# solve(K)
# a <- get_a_vec_3(l0,m0,m,r)
# w <- get_W_3(a,r,m0)
# sum(w)-(r+1)
# for(k in 1:l){
#   temp <- 0
#   for(i in (-r):0){
#     temp <- temp + (i^k)*(w[i+r+1])
#     #print(temp)
#   }
#   print(temp)
# }


#
# a <- get_coefs(l0,m0,m,r)
# w <- get_W(a,r,m0)
# sum(w)-(2*r+1) # this should be 2r+1
# # check the second condition also
# for(k in 1:l){
#   temp <- 0
#   for(i in -r:r){
#     temp <- temp + (i^k)*(w[-i+r+1])
#   }
#   print(temp)
# }

# these seem to be correct for all r now
