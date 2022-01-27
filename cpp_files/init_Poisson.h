#include "onlineblock.h"
// #include "blockPoisson.h"


// [[Rcpp::export]]
arma::mat update_init_Tau(arma::mat tau, int m, int K){
  arma::mat tau_mod;
  tau_mod.zeros(m,K);
  for(int i =0; i < m; ++i){
    tau_mod.row(i) = correct_tau(tau.row(i));
  }
  return tau_mod;
}

// full online estimation procedure with initial values
// [[Rcpp::export]]
Rcpp::List estimate_Poisson_init(
    arma::mat full_data, 
    Rcpp::List A, //arma::mat A,
    int m,
    int K,
    double T,
    double dT,
    arma::mat B,
    int inter_T,
    arma::mat tau, // pass this in, from the init procedure
    double start,
    bool is_elbo = false
){
  // iterate this over the time windows...
  int N = int((T - start)/dT); // account for init
  int slices = int(N/inter_T);
  double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  int ind = 0;
  int nall = full_data.n_rows;
  arma::rowvec Pi(K);
  Pi.fill(1.0 / K);
  arma::mat S(m,K);
  // arma::mat tau(m,K);
  S.fill(1.0/K);
  // tau.fill(1.0/K);
  tau = update_init_Tau(tau, m, K);
  arma::cube inter_tau(m,K,slices+1);
  arma::cube inter_B(K,K,N);
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll;
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
  ave_elbo.zeros(N);
  int cum_events = 0;
  cout<<"Gotten to here \n"<<endl;
  for(int n = 0; n < N; ++n){
    cout<<n<<endl;
    double Tn = dT*(n+1) + start; //account for init
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
    eta = 1/pow(1 + n, .5)/sub_data.n_rows*(K*K);
    S = updateS(sub_data, tau, B, A, S, K, m, dT);
    tau = updateTau(S,Pi,m,K); 
    B = updateB(sub_data, tau, B, K, A, m, dT, eta);
    inter_B.slice(n) = B;
    Pi = updatePi(tau,K);
    if (is_elbo) {
      // changed these to just the current window data
      curr_elbo(n) = computeELBO(sub_data, tau, B, Pi, A, m, K, dT);
      ave_elbo(n) = curr_elbo(n)/cum_events;
      curr_ll(n) = computeLL(sub_data, tau, B, Pi, A, m, K, t_curr);
      ave_ll(n) = curr_ll(n)/cum_events;
    }
    if(n % inter_T == 0 ){
      inter_tau.slice(ind) = tau;
      ind = ind + 1;
    }
    
  }
  
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            Named("early_tau")= inter_tau,
                            Named("inter_B") = inter_B,
                            Named("B")=B,
                            Named("Pi")=Pi,
                            Named("wind_elbo") = curr_elbo,
                            Named("AveELBO")=ave_elbo,
                            Named("wind_ll") = curr_ll,
                            Named("logL") = ave_ll);
}
