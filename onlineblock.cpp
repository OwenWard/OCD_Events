#include "onlineblock.h"
#include "supportsim.h"




#include "blockhawkes.h"
#include "blockPoisson.h"

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

