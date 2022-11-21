#include "onlineblock.h"
// general function to compute log likelihood for inhomogeneous processes


Rcpp::List compute_regret_inhom(
    arma::mat full_data,
    Rcpp::List A, //arma::mat A,
    int m,
    int K,
    int H,
    double T,
    double dT,
    double window,
    arma::vec true_z,
    arma::mat MuA_ests,
    // will modify this, each row will be a flattened cube
    arma::cube tau_ests, // will this actually be used?
    arma::cube true_MuA // and this, to cube
){
  // convert the true z vector to a tau matrix,
  // then can use compute llhood function
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
  //
  int cum_events = 0;
  double t_prev, t_length;
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
    elbo_dat = full_data.rows(0,end_pos);
    start_pos = curr_pos;
    t_length = t_curr-t_prev;
    if(n > 0) {
      arma::rowvec curr = MuA_ests.row(n-1);
      arma::cube prev_MuA(K, K, H);
      prev_MuA = to_cube(curr, K, K, H);
      arma::mat prev_tau = tau_ests.slice(n-1);
      // TO DO: make these tau's exact, i.e 0 or 1
      arma::rowvec prev_Pi;
      prev_Pi = sum(prev_tau, 0)/m;

      pred_ll(n-1) = get_elbo_nonhomoHak(sub_data, t_prev,
                                t_curr, prev_tau, prev_MuA, B,
                                prev_Pi, A, lam, m, K, H, window);
      
      ave_pred_ll(n-1) = pred_ll(n-1)/sub_data.n_rows;
      // compute same for batch estimator...
      
    }
    arma::rowvec curr = MuA_ests.row(n);
    arma::cube curr_MuA(K, K, H);
    curr_MuA = to_cube(curr, K, K, H);
    arma::mat curr_tau = tau_ests.slice(n);
    // TO DO, make these exact also
    arma::rowvec curr_Pi;
    curr_Pi = sum(curr_tau, 0)/m;
    // likelihood using known tau
    cout<<sub_data.n_rows<<endl;
    // TO DO, check tau for here
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
  
}


// // [[Rcpp::export]]
// Rcpp::List compute_regret_inhom(
//     arma::mat full_data,
//     Rcpp::List A, //arma::mat A,
//     int m,
//     int K,
//     int H,
//     double T,
//     double dT,
//     double window,
//     arma::vec true_z,
//     arma::cube MuA_ests,
//     // will need to modify this, will now be seq of 3d arrays
//     arma::cube tau_ests,
//     arma::mat true_B // and this, to cube
// ){
//   // convert the true z vector to a tau matrix,
//   // then can use compute llhood function
//   // arma::mat tau(m,K);
//   // tau.fill(0);
//   // // then iterate over z
//   // for(int i=0; i<m; ++i){
//   //   int ind = true_z[i]-1;
//   //   tau(i,ind) = 1;
//   // }
//   // arma::rowvec Pi;
//   // Pi = sum(tau, 0)/m;
//   // // Pi.print("Pi:");
//   // int N = int(T/dT);
//   // // double eta;
//   // int start_pos = 0;
//   // int curr_pos = 0;
//   // int end_pos = 0;
//   // // int ind = 0;
//   // int nall = full_data.n_rows;
//   // arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, true_ll, av_true_ll, pred_ll;
//   // arma::vec event_loss_trueZ, event_loss_estZ, ave_pred_ll, emp_ll;
//   // event_loss_trueZ.zeros(N);
//   // event_loss_estZ.zeros(N);
//   // curr_elbo.zeros(N);
//   // curr_ll.zeros(N);
//   // ave_ll.zeros(N);
//   // true_ll.zeros(N);
//   // emp_ll.zeros(N);
//   // av_true_ll.zeros(N);
//   // ave_elbo.zeros(N);
//   // pred_ll.zeros(N-1);
//   // ave_pred_ll.zeros(N-1);
// 
//   // for Poisson B is fixed to be zero
//   // arma::mat B(K, K);
//   // B.fill(0.0);
//   // //
//   // int cum_events = 0;
//   // double t_prev, t_length;
//   // double lam = 0.0;
//   // t_prev = 0;
//   // for(int n = 0; n < N; ++n){
//   //   double Tn = dT*(n+1);
//   //   arma::rowvec event = full_data.row(start_pos);
//   //   double t_curr = event(2);
//   //   while(t_curr <= Tn){
//   //     if(curr_pos >= nall-1){
//   //       break;
//   //     }
//   //     else{
//   //       curr_pos += 1;
//   //       event = full_data.row(curr_pos);
//   //       t_curr = event(2);
//   //     }
//   //   }
//   //   end_pos = curr_pos;
//   //   arma::mat sub_data, elbo_dat;
//   //   sub_data = full_data.rows(start_pos, end_pos);
//   //   cum_events += sub_data.n_rows;
//   //   elbo_dat = full_data.rows(0,end_pos);
//   //   // this is the cumulative data, should it be just the events in that
//   //   // window here?
//   //   start_pos = curr_pos;
//   //   t_length = t_curr-t_prev;
//   //   // compute predictive log likelihood here
//   //   if(n > 0) {
//   //     // TO DO, update this
//   //     // arma::mat prev_MuA = MuA_ests.slice(n-1);
//   //     // // reshape this here
//   //     // arma::mat prev_tau = tau_ests.slice(n-1);
//   //     // arma::rowvec prev_Pi;
//   //     // prev_Pi = sum(prev_tau, 0)/m;
//   // 
//   //     // need to estimate z_prev, make it R indexed, or does ll_nonhomo do this?
//   //     // define B to be 0 here somewhere
//   // 
//   //     //
//   //     // unordered_map<string, arma::vec> datamap = transfer(sub_data);
//   // 
//   // 
//   //     // pred_ll(n-1) = ll_nonhomo(datamap, t_prev,
//   //     //                           t_curr, prev_tau, prev_MuA, B,
//   //     //                           prev_Pi, A, lam, m, K, H, window);
//   //       // (sub_data, prev_tau, prev_B, prev_Pi,
//   //       //       A, m, K, t_length);
//   //     // ave_pred_ll(n-1) = pred_ll(n-1)/sub_data.n_rows;
//   //     // compute same for batch estimator...
//   // 
//   //   }
//   //   t_prev = t_curr;
//   //   arma::mat curr_MuA = MuA_ests.slice(n);
//   //   // reshape this here
//   //   arma::mat curr_tau = tau_ests.slice(n);
//   //   arma::rowvec curr_Pi;
//   //   curr_Pi = sum(curr_tau, 0)/m;
//   //   // likelihood using known tau
//     // cout<<sub_data.n_rows<<endl;
// 
//     // curr_ll(n) = ll_nonhomo(sub_data, t_prev,
//     //                         t_curr, curr_tau, curr_MuA, B,
//     //                         curr_Pi, A, lam, m, K, H, window);
//     // curr_ll(n) = computeLL(sub_data, tau, curr_B, curr_Pi, A, m, K, t_length);
// 
//     // ave_ll(n) = curr_ll(n)/cum_events;
// 
//     // TO DO, convert true_z to corresponding tau representation
//     // arma::mat tau_true(m, K);
//     // tau_true.fill(0.0);
//     // arma::mat MuA_true(K, K);
//     // MuA_true.fill(0.0);
// 
//     // true_ll(n) = ll_nonhomo(sub_data, t_prev, t_curr, tau_true, MuA_true, B,
//     //                         curr_Pi, A, lam, m, K, H, window);
//     // // also compute "empirical regret" here,
//     // // using estimated z and theta
//     // emp_ll(n) = ll_nonhomo(sub_data, t_prev, t_curr, curr_tau, curr_MuA,
//     //                        B, curr_Pi, A, lam, m, K, H, window);
//     //
// 
//     // was sub_data
//     // av_true_ll(n) = true_ll(n)/cum_events;
//     //
// 
//     // TO DO, update the below after writing the functions for them
// 
//     // event_loss_trueZ(n) = curr_loss(sub_data,
//     //                  curr_B,
//     //                  true_z,
//     //                  curr_tau, K, dT, m)["True_z"];
//     // // need to pass inter tau in here
//     // event_loss_estZ(n) = curr_loss(sub_data,
//     //                 curr_B, true_z, curr_tau, K, dT, m)["Est_z"];
//   // }
//   // return Rcpp::List::create(Named("EstLLH") = curr_ll,
//   //                           Named("Ave_est_LLH") = ave_ll,
//   //                           Named("TrueLLH") = true_ll,
//   //                           Named("Ave_true_LLH") = av_true_ll,
//   //                           Named("EmpLLH") = emp_ll,
//   //                           Named("Pred_LL") = pred_ll,
//   //                           Named("Ave_Pred_LL") = ave_pred_ll,
//   //                           Named("Online_Loss_True") = event_loss_trueZ,
//   //                           Named("Online_Loss_Est") = event_loss_estZ);
//   return Rcpp::List::create(Named("Test1") = 0.01,
//                             Named("Test2") = 0.02);
// }


