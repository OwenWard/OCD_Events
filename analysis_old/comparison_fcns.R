library(mclust)
library(dplyr)
library(ppsbm)
library(kernlab)
library(blockmodels)
library(Rcpp)
library(RcppArmadillo)
# sourceCpp("onlineblock.cpp")



spec_ARI <- function(a_matrix,K, true_Z){
  sc <- specc(a_matrix, centers=K)
  estimated_groups = sc@.Data
  return(adjustedRandIndex(estimated_groups,true_Z))
}


sbm_bin_ARI <- function(bin_matrix,K, true_Z){
  sbm_bin = BM_bernoulli('SBM',bin_matrix,explore_min = K, explore_max = K,
                         verbosity = 0)
  sbm_bin$estimate()
  sbm_bin_est = apply(sbm_bin$memberships[[3]]$Z,1,which.max)
  return(adjustedRandIndex(true_Z,sbm_bin_est))
}


sbm_count_ARI <- function(count_matrix,K,true_Z){
  sbm_count = BM_poisson('SBM',count_matrix, explore_min = K, explore_max = K,
                         verbosity = 0)
  sbm_count$estimate()
  sbm_count_est = apply(sbm_count$memberships[[3]]$Z,1,which.max)
  return(adjustedRandIndex(true_Z,sbm_count_est))
}
