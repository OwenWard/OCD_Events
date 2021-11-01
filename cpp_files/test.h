#include "onlineblock.h"

// for testing out some c++ code that may be added in as needed

// [[Rcpp::export]]
Rcpp::List init_vals(int K, int n){
  arma::rowvec Pi(K);
  Pi.fill(1.0/K);
  arma::mat tau(n,K);
  tau.fill(1.0/K);
  arma::mat S(n,K);
  S.fill(1.0/K);
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            Named("Pi")=Pi);
}



// [[Rcpp::export]]
Rcpp::List compute_regret(
    arma::mat full_data, 
    Rcpp::List A, //arma::mat A,
    int m,
    int K,
    double T,
    double dT,
    arma::vec true_z,
    arma::cube B_ests,
    arma::mat true_B
){
  // convert the true z vector to a tau matrix,
  // then can use compute llhood function
  arma::mat tau(m,K);
  tau.fill(0);
  // then iterate over z
  for(int i=0; i<m; ++i){
    int ind = true_z[i]-1;
    tau(i,ind) = 1;
  }
  arma::rowvec Pi;
  Pi = sum(tau, 0)/m;
  // Pi.print("Pi:");
  int N = int(T/dT);
  // double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  // int ind = 0;
  int nall = full_data.n_rows;
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, true_ll, av_true_ll;
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
  true_ll.zeros(N);
  av_true_ll.zeros(N);
  ave_elbo.zeros(N);
  int cum_events = 0;
  for(int n = 0; n < N; ++n){
    double Tn = dT*(n+1);
    arma::rowvec event = full_data.row(start_pos);
    double t_curr = event(2);
    while(t_curr <= Tn){
      if(curr_pos >= nall-1){
        break;
      }
      else{
        curr_pos += 1;
        event = full_data.row(curr_pos);
        t_curr = event(2);
      }
    }
    end_pos = curr_pos;
    arma::mat sub_data, elbo_dat;
    sub_data = full_data.rows(start_pos, end_pos);
    cum_events += sub_data.n_rows;
    elbo_dat = full_data.rows(0,end_pos); 
    // this is the cumulative data, should it be just the events in that
    // window here?
    start_pos = curr_pos;
    // cout<<"here"<<endl;
    // R.printf();
    // Rcpp::print(B_ests.slice(n));
    arma::mat curr_B = B_ests.slice(n);
    curr_ll(n) = computeLL(sub_data, tau, curr_B, Pi, A, m, K, t_curr);
    ave_ll(n) = curr_ll(n)/cum_events;
    true_ll(n) = computeLL(sub_data, tau, true_B, Pi, A, m, K, t_curr);
    av_true_ll(n) = true_ll(n)/cum_events;
  }
  return Rcpp::List::create(Named("EstLLH") = curr_ll,
                            Named("Ave_est_LLH") = ave_ll,
                            Named("TrueLLH") = true_ll,
                            Named("Ave_true_LLH") = av_true_ll);
}

