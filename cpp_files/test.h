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
Rcpp::List curr_loss(
    arma::mat curr_data,
    arma::mat curr_B,
    arma::vec true_z, // computed in R, so not zero indexed
    arma::mat curr_tau, 
    int K,
    double dT,
    int m) {
  // here we compute the loss for the number of events in a window
  // of length dT, given the estimates for B and tau after each of these
  // windows (already computed)
  //
  
  int n = curr_data.n_rows;
  arma::mat event_counts(m,m), diff_mat(m,m), diff_mat2(m,m);
  diff_mat.fill(0); // difference given true z
  diff_mat2.fill(0); // difference given estimated z

  for(int i = 0; i < n; ++i) {
    arma::rowvec event = curr_data.row(i);
    int send = event(0);
    int rec = event(1);
    event_counts(send, rec) += 1;
  }

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
  norm_diff = abs(norm_diff);
  // this is the loss with the true z estimates
  
  // estimate z given tau here, repeat for those
  arma::vec est_z;
  est_z.zeros(m);
  for(int i = 0; i < m; ++i){
    arma::rowvec tau_i = curr_tau.row(i);
    est_z(i) = tau_i.index_max();
  }
  for(int i = 0; i < m; ++i){
    for(int j = 0; j < m; ++j){
      if(i != j ){
        int send_ind = est_z(i); // for zero index
        int rec_ind = est_z(j);
        double curr_rate = dT * curr_B(send_ind, rec_ind);
        diff_mat2(i,j) = event_counts(i,j) - curr_rate;
      }
    }
  }
  // without knowning the true z
  double norm_diff2 = accu(diff_mat2)/ pow(m, 2);
  norm_diff2 = abs(norm_diff2);
  return Rcpp::List::create(Named("True_z") = norm_diff,
                            Named("Est_z") = norm_diff2);
}


// [[Rcpp::export]]
Rcpp::List batch_loss(
  arma::mat full_data,
  Rcpp::List A,
  int m,
  int K,
  double T,
  double dT,
  arma::mat batch_B,
  arma::vec true_z
){
  // compute the loss for each window
  // using the estimated z and batch B each time
  // along with batch predicted log likelihood also
  int N = int(T/dT);
  // double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  // int ind = 0;
  int nall = full_data.n_rows;
  
  arma::vec batch_loss, batch_pred_ll, batch_ave_pred_ll;
  arma::mat tau(m,K); //used in pred llh
  batch_loss.zeros(N);
  batch_pred_ll.zeros(N-1);
  batch_ave_pred_ll.zeros(N-1);
  tau.fill(0);
  // then iterate over z
  for(int i=0; i<m; ++i){
    int ind = true_z[i]-1;
    tau(i,ind) = 1;
  }
  arma::rowvec Pi;
  Pi = sum(tau, 0)/m;
  // need to compute the length of the window each time here...
  double t_prev, t_length;
  t_prev = 0;
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
    start_pos = curr_pos;
    // then compute the loss here
    batch_loss(n) = curr_loss(sub_data,
              batch_B,
              true_z,
              tau,
              K, dT, m)["True_z"];
    if(n > 0) {
      // need to actually get the tau for this I guess
      t_length = t_curr - t_prev;
      batch_pred_ll(n-1) = computeLL(sub_data, tau, batch_B, Pi, A,
                    m, K, t_length);
      batch_ave_pred_ll(n-1) = batch_pred_ll(n-1)/sub_data.n_rows;
      // cout<<sub_data.n_rows<<endl;
    }
    t_prev = t_curr;
  }
  return Rcpp::List::create(Named("Batch_loss") = batch_loss,
                            Named("Batch_Pred_LL") = batch_pred_ll,
                            Named("Batch_Ave_Pred_LL") = batch_ave_pred_ll);
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
    arma::cube tau_ests,
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
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, true_ll, av_true_ll, pred_ll;
  arma::vec event_loss_trueZ, event_loss_estZ, ave_pred_ll;
  event_loss_trueZ.zeros(N);
  event_loss_estZ.zeros(N);
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
  true_ll.zeros(N);
  av_true_ll.zeros(N);
  ave_elbo.zeros(N);
  pred_ll.zeros(N-1);
  ave_pred_ll.zeros(N-1);
  int cum_events = 0;
  double t_prev, t_length;
  t_prev = 0;
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
    
    // compute predictive log likelihood here
    if(n > 0) {
      arma::mat prev_B = B_ests.slice(n-1);
      arma::mat prev_tau = tau_ests.slice(n-1);
      arma::rowvec prev_Pi;
      prev_Pi = sum(prev_tau, 0)/m;
      t_length = t_curr-t_prev;
      pred_ll(n-1) = computeLL(sub_data, prev_tau, prev_B, prev_Pi,
              A, m, K, t_length);
      ave_pred_ll(n-1) = pred_ll(n-1)/sub_data.n_rows;
      // compute same for batch estimator...
      
    }
    t_prev = t_curr;
    
    ///
    arma::mat curr_B = B_ests.slice(n);
    arma::mat curr_tau = tau_ests.slice(n);
    arma::rowvec curr_Pi;
    curr_Pi = sum(curr_tau, 0)/m;
    // likelihood using known tau
    curr_ll(n) = computeLL(sub_dat, tau, curr_B, curr_Pi, A, m, K, t_curr);
    ave_ll(n) = curr_ll(n)/cum_events;
    true_ll(n) = computeLL(sub_dat, tau, true_B, curr_Pi, A, m, K, t_curr);
    // was sub_data
    av_true_ll(n) = true_ll(n)/cum_events;
    //
    event_loss_trueZ(n) = curr_loss(sub_data,
                     curr_B,
                     true_z,
                     curr_tau, K, dT, m)["True_z"];
    // need to pass inter tau in here
    event_loss_estZ(n) = curr_loss(sub_data,
                    curr_B, true_z, curr_tau, K, dT, m)["Est_z"];
  }
  return Rcpp::List::create(Named("EstLLH") = curr_ll,
                            Named("Ave_est_LLH") = ave_ll,
                            Named("TrueLLH") = true_ll,
                            Named("Ave_true_LLH") = av_true_ll,
                            Named("Pred_LL") = pred_ll,
                            Named("Ave_Pred_LL") = ave_pred_ll,
                            Named("Online_Loss_True") = event_loss_trueZ,
                            Named("Online_Loss_Est") = event_loss_estZ);
}
