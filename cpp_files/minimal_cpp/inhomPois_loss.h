#include "onlineblock.h"
// general function to compute log likelihood for inhomogeneous processes


// [[Rcpp::export]]
Rcpp::List inhom_batch_loss(
    arma::mat full_data,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    arma::cube batch_B,
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
Rcpp::List inhom_curr_loss(
    arma::mat curr_data,
    arma::cube curr_MuA,
    arma::vec true_z, // computed in R, so not zero indexed
    arma::mat curr_tau, 
    int K,
    double dT,
    double t_start,
    double t_end, 
    double window,
    int H,
    int m) {
  // here we compute the loss for the number of events in a window
  // of length dT, given the estimates for B and tau after each of these
  // windows (already computed)
  //
  
  arma::vec tvec(H);
  tvec.fill(0.0);
  int h1 = floor(t_start/window);
  int h2 = floor(t_end/window);
  
  for (int w = h1; w <= h2; w++){
    int h = w % H;
    if (h1 == h2){
      tvec(h) = t_end - t_start;
      break;
    }
    if (w == h1) {
      tvec(h) += (w + 1)*window - t_start;
    } else if (w == h2) {
      tvec(h) += t_end - w * window;
    } else {
      tvec(h) += window;
    }
  }
  
  
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
        double curr_rate = 0.0;
        for (int h = 0; h < H; h++) {
          curr_rate += curr_MuA(send_ind, rec_ind, h) * tvec(h);
        }
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
        double curr_rate = 0.0;
        for(int h = 0; h < H; h++) {
          curr_rate += curr_MuA(send_ind, rec_ind, h) * tvec(h);
        }
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
Rcpp::List compute_regret_inhom(arma::mat full_data,
    Rcpp::List A, //arma::mat A,
    int m,
    int K,
    int H,
    double T,
    double dT,
    double window,
    arma::vec true_z,
    arma::mat MuA_ests, // will modify this, each row will be a flattened cube
    arma::cube tau_ests, // will this actually be used?
    arma::cube true_MuA){ // and this, to cube
  arma::mat tau_true(m,K);
  tau_true.fill(0);
  // then iterate over z
  for(int i=0; i<m; ++i){
    int ind = true_z[i]-1;
    tau_true(i,ind) = 1;
  }
  arma::rowvec Pi_true;
  Pi_true = sum(tau_true, 0)/m;

  int N = int(T/dT);
  // double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  // int ind = 0;
  int nall = full_data.n_rows;
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, true_ll, av_true_ll, pred_ll;
  arma::vec event_loss_trueZ, event_loss_estZ, ave_pred_ll, emp_ll;
  event_loss_trueZ.zeros(N);
  event_loss_estZ.zeros(N);
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
  true_ll.zeros(N);
  emp_ll.zeros(N);
  av_true_ll.zeros(N);
  ave_elbo.zeros(N);
  pred_ll.zeros(N-1);
  ave_pred_ll.zeros(N-1);
  // for Poisson B is fixed to be zero
  arma::mat B(K, K);
  B.fill(0.0);
  int cum_events = 0;
  double t_prev;
  double lam = 0.0;
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
    
    // elbo_dat = full_data.rows(0,end_pos);
    start_pos = curr_pos;
    if(n > 0) {
      arma::rowvec curr = MuA_ests.row(n-1);
      arma::cube prev_MuA(K, K, H);
      prev_MuA = to_cube(curr, K, K, H);
      arma::mat prev_tau(m, K);
      prev_tau = tau_ests.slice(n-1);
      prev_tau = max_tau(prev_tau);
      //make these tau's exact, i.e 0 or 1
      arma::rowvec prev_Pi;
      prev_Pi = sum(prev_tau, 0)/m;

      pred_ll(n-1) = get_elbo_nonhomoHak(sub_data, t_prev,
                                t_curr, prev_tau, prev_MuA, B,
                                prev_Pi, A, lam, m, K, H, window);

      ave_pred_ll(n-1) = pred_ll(n-1)/sub_data.n_rows;
      // compute same for batch estimator...

    }
    cout<<n<<endl;
    // this is causing the problem here
    arma::rowvec curr = MuA_ests.row(n);
    arma::cube curr_MuA(K, K, H);
    curr_MuA = to_cube(curr, K, K, H);
    arma::mat curr_tau = tau_ests.slice(n);
    curr_tau = max_tau(curr_tau);
    arma::rowvec curr_Pi;
    curr_Pi = sum(curr_tau, 0)/m;
    // likelihood using known tau
    cout<<sub_data.n_rows<<endl;
    curr_ll(n) = get_elbo_nonhomoHak(sub_data, t_prev,
                            t_curr, curr_tau, curr_MuA, B,
                            curr_Pi, A, lam, m, K, H, window);
    ave_ll(n) = curr_ll(n)/cum_events;


    true_ll(n) = get_elbo_nonhomoHak(sub_data, t_prev, t_curr, tau_true,
                                     true_MuA, B,
                                     Pi_true, A, lam, m, K, H, window);
    // also compute "empirical regret" here,
    // // using estimated z and theta
    emp_ll(n) = get_elbo_nonhomoHak(sub_data, t_prev, t_curr, curr_tau,
                                    curr_MuA,
                                    B, curr_Pi, A, lam, m, K, H, window);
    av_true_ll(n) = true_ll(n)/cum_events;
    
    event_loss_trueZ(n) = inhom_curr_loss(sub_data,
                     curr_MuA,
                     true_z,
                     curr_tau, K, dT, t_prev,
                     t_curr, window, H, m)["True_z"];
    // need to pass inter tau in here
    event_loss_estZ(n) = inhom_curr_loss(sub_data,
                     curr_MuA, true_z, curr_tau,
                     K, dT, t_prev, t_curr, window, H, m)["Est_z"];
    
    t_prev = t_curr;
  }

  return Rcpp::List::create(Named("EstLLH") = curr_ll,
                            Named("Ave_est_LLH") = ave_ll,
                            Named("TrueLLH") = true_ll,
                            Named("Ave_true_LLH") = av_true_ll,
                            Named("EmpLLH") = emp_ll,
                            Named("Pred_LL") = pred_ll,
                            Named("Ave_Pred_LL") = ave_pred_ll,
                            Named("Online_Loss_True") = event_loss_trueZ,
                            Named("Online_Loss_Est") = event_loss_estZ);
  // return 0;
}
