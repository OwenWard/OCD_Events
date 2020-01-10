#include "onlineblock.h"


// [[Rcpp::export]]
arma::mat Predict_Counts_Hawkes(
    double startT,
    double finalT,
    Rcpp::List A,
    arma::vec Z,
    arma::mat Mu,
    arma::mat B,
    int m,
    double lam){
  arma::mat counts;
  int z1, z2;
  double mu,b,exp_term,k;
  counts.zeros(1,3);
  double mean_events;
  arma::rowvec curr_counts;
  curr_counts.zeros(3);
  int n_edge;
  for(int i = 0; i < m; i++){
    arma::rowvec edge = A[i];
    n_edge = edge.n_elem;
    for(int p = 0; p < n_edge; p++){
    //for(int j = 0; j < m; j++){
      int j = edge(p);
      if(i != j){
        z1 = Z(i), z2 = Z(j);
        mu = Mu(z1, z2), b = B(z1, z2);
        k = lam*(1-b);
        exp_term = exp(-k*startT) - exp(-k*finalT);
        mean_events = (finalT-startT)-b/(lam*pow(1-b,2))*exp_term;
        mean_events = mean_events*mu/(1-b);
        // create a vector storing i j and mean events
        int nrows = counts.n_rows;
        curr_counts(0) = i;
        curr_counts(1) = j;
        curr_counts(2) = mean_events;
        counts.insert_rows(nrows,curr_counts);
      }
    }
  }
  return counts;
}


