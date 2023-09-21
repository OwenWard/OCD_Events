#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "onlineblock.h"
#include "supportsim.h"



// #include "blockPoisson.h"
// #include "blockhawkes.h"
// #include "nonhomohak.h"
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

// // [[Rcpp::export]]
// Rcpp::List test(
// 	arma::mat alltimes,
// 	arma::mat tau,
// 	arma::mat Mu,
// 	arma::mat B,
// 	arma::rowvec Pi,
// 	arma::mat S,
// 	double t_start,
// 	double Tn,
// 	int m,
// 	int K,
// 	Rcpp::List A,
// 	double lam,
// 	double eta
// 	){
// 	unordered_map<string, arma::vec> datamap = transfer(alltimes);
// 	Rcpp::List paralist = update_on(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
// 	return paralist;
// }
// 

