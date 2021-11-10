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
double curr_loss(
    arma::mat curr_data,
    arma::mat curr_B,
    arma::vec true_z,
    int K,
    double dT,
    int m) {
  int n = curr_data.n_rows;
  arma::mat event_counts(m,m), diff_mat(m,m);
  // double total = 0;
  
  for(int i = 0; i < n; ++i) {
    arma::rowvec event = curr_data.row(i);
    int send = event(0);
    int rec = event(1);
    event_counts(send, rec) += 1;
  }
  // return event_counts;
  // this is correct
  // 
  // cout<<"Here"<<endl;
  for(int i = 0; i < m; ++i){
    for(int j = 0; j < m; ++j){
      if(i != j ){
        int send_ind = true_z(i) - 1; // for zero index
        int rec_ind = true_z(j) - 1;
        double curr_rate = dT * curr_B(send_ind, rec_ind);
        diff_mat(i,j) = event_counts(i,j) - curr_rate;
      }
    }
  }
  double norm_diff = accu(diff_mat)/ pow(m, 2); 
  return abs(norm_diff).s;
  // want to sum this and then divide by n^2?
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
  arma::vec event_loss;
  event_loss.zeros(N);
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
    curr_ll(n) = computeLL(elbo_dat, tau, curr_B, Pi, A, m, K, t_curr);
    ave_ll(n) = curr_ll(n)/cum_events;
    true_ll(n) = computeLL(elbo_dat, tau, true_B, Pi, A, m, K, t_curr);
    // was sub_data
    av_true_ll(n) = true_ll(n)/cum_events;
    //
    event_loss(n) = curr_loss(sub_data, curr_B, true_z, K, dT, m);
    
  }
  return Rcpp::List::create(Named("EstLLH") = curr_ll,
                            Named("Ave_est_LLH") = ave_ll,
                            Named("TrueLLH") = true_ll,
                            Named("Ave_true_LLH") = av_true_ll,
                            Named("Online_Loss")= event_loss);
}



// // [[Rcpp::export]]
// double curr_loss(
//   arma::mat curr_data,
//   // arma::mat tau,/
//   // arma::mat curr_B,
//   // double dT,
//   int K,
//   int m
// ){
//   // construct event matrix and expected rate matrix
//   arma::mat event_counts(m,m);
//   event_counts.fill(0);
//   // int num_events = curr_data.n_rows;
//   // for(int i = 0; i < num_events; ++i){
//   //   arma::rowvec event = curr_data.row(i);
//   //   int send = event(0);
//   //   int rec = event(1);
//   //   event_counts(send, rec) += 1;
//   // }
//   //event_counts.print();
//   //return Rcpp::List::create(Named("events") = event_counts);
  // double diff = 1.0;
  // return 1.0;
// }


// loss in terms of predicted events over time...
// assuming known true z vectors say!
// Rcpp::List compute_loss(
//     arma::mat full_data, 
//     Rcpp::List A, //arma::mat A,
//     arma::mat tau,
//     int m,
//     int K,
//     double T,
//     double dT,
//     arma::vec true_z,
//     arma::cube B_ests,
// ){
//   // convert the true z vector to a tau matrix,
//   // then can use compute llhood function
//   arma::mat tau(m,K);
//   tau.fill(0);
//   // then iterate over z
//   for(int i=0; i<m; ++i){
//     int ind = true_z[i]-1;
//     tau(i,ind) = 1;
//   }
//   arma::rowvec Pi;
//   Pi = sum(tau, 0)/m;
//   // Pi.print("Pi:");
//   int N = int(T/dT);
//   // double eta;
//   int start_pos = 0;
//   int curr_pos = 0;
//   int end_pos = 0;
//   // int ind = 0;
//   int nall = full_data.n_rows;
//   arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, true_ll, av_true_ll;
//   curr_elbo.zeros(N);
//   curr_ll.zeros(N);
//   ave_ll.zeros(N);
//   true_ll.zeros(N);
//   av_true_ll.zeros(N);
//   ave_elbo.zeros(N);
//   int cum_events = 0;
//   for(int n = 0; n < N; ++n){
//     double Tn = dT*(n+1);
//     arma::rowvec event = full_data.row(start_pos);
//     double t_curr = event(2);
//     while(t_curr <= Tn){
//       if(curr_pos >= nall-1){
//         break;
//       }
//       else{
//         curr_pos += 1;
//         event = full_data.row(curr_pos);
//         t_curr = event(2);
//       }
//     }
//     end_pos = curr_pos;
//     arma::mat sub_data, elbo_dat;
//     sub_data = full_data.rows(start_pos, end_pos);
//     // then want to compare predicted events vs true, using sub_data
//     
//     cum_events += sub_data.n_rows;
//     elbo_dat = full_data.rows(0,end_pos); 
//     // this is the cumulative data, should it be just the events in that
//     // window here?
//     start_pos = curr_pos;
//     // cout<<"here"<<endl;
//     // R.printf();
//     // Rcpp::print(B_ests.slice(n));
//     arma::mat curr_B = B_ests.slice(n);
//     curr_ll(n) = computeLL(elbo_dat, tau, curr_B, Pi, A, m, K, t_curr);
//     ave_ll(n) = curr_ll(n)/cum_events;
//     true_ll(n) = computeLL(elbo_dat, tau, true_B, Pi, A, m, K, t_curr);
//     // was sub_data
//     av_true_ll(n) = true_ll(n)/cum_events;
//   }
//   return Rcpp::List::create(Named("EstLLH") = curr_ll,
//                             Named("Ave_est_LLH") = ave_ll,
//                             Named("TrueLLH") = true_ll,
//                             Named("Ave_true_LLH") = av_true_ll);
// }
