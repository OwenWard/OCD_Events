#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "onlineblock.h"
#include "supportsim.h"


#include "blockPoisson.h"
#include "blockhawkes.h"
#include "inhomPoisson.h" 
#include "inhomHawkes.h"
// #include "link_predict.h"
// #include "ccrm.h"
// #include "new_predict.h"
// #include "test.h"
// #include "init_Poisson.h"

// [[Rcpp::export]]
arma::vec func(arma::vec X){
  X(0) = 4;
  arma::vec y;
  return y;
}

// [[Rcpp::export]]
void test(){
  arma::vec X(3);
  X.fill(0.0);
  X = func(X);
  X.print();
}



