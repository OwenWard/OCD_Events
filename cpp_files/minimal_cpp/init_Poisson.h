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
  arma::cube inter_tau(m,K,slices);
  arma::cube inter_B(K,K,N);
  arma::cube inter_S(m, K, N);
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, full_elbo;
  arma::vec cum_elbo, cum_events_vec;
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
  ave_elbo.zeros(N);
  full_elbo.zeros(N);
  cum_elbo.zeros(N);
  cum_events_vec.zeros(N);
  int cum_events = 0;
  for(int n = 0; n < N; ++n){
    // cout<<n<<endl;
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
    // cout << sub_data.n_rows << endl;
    cum_events += sub_data.n_rows;
    elbo_dat = full_data.rows(0,end_pos); 
    start_pos = curr_pos;
    eta = 1/pow(1 + n, .5)/sub_data.n_rows*(K*K);
    S = updateS(sub_data, tau, B, A, S, K, m, dT);
    // cout << S << endl;
    inter_S.slice(n) = S;
    tau = updateTau(S,Pi,m,K); 
    B = updateB(sub_data, tau, B, K, A, m, dT, eta);
    inter_B.slice(n) = B;
    Pi = updatePi(tau,K);
    // if(n == 0){
    //   cout << S << endl;
    // }
    if (is_elbo) {
      // changed these to just the current window data
      curr_elbo(n) = computeELBO(sub_data, tau, B, Pi, A, m, K, dT);
      ave_elbo(n) = curr_elbo(n)/cum_events;
      curr_ll(n) = computeLL(sub_data, tau, B, Pi, A, m, K, t_curr);
      ave_ll(n) = curr_ll(n)/cum_events;
      full_elbo(n) = computeELBO(full_data, tau, B, Pi, A, m, K, T);
      cum_elbo(n) = computeELBO(elbo_dat, tau, B, Pi, A, m, K, t_curr)/cum_events;
      cum_events_vec(n) = cum_events;
    }
    if(n % inter_T == 0 ){
      inter_tau.slice(ind) = tau;
      ind = ind + 1;
    }
    
  }
  
  return Rcpp::List::create(Named("S")= inter_S,
                            Named("tau")=tau,
                            Named("early_tau")= inter_tau,
                            Named("inter_B") = inter_B,
                            Named("B")=B,
                            Named("Pi")=Pi,
                            Named("wind_elbo") = curr_elbo,
                            Named("AveELBO")=ave_elbo,
                            Named("wind_ll") = curr_ll,
                            Named("logL") = ave_ll,
                            Named("full_ELBO") = full_elbo,
                            Named("Cum_ELBO") = cum_elbo,
                            Named("Cum_Events") = cum_events_vec);
}





//// modify inhomogeneous Poisson to take initial values also
//// pass in the community assignments and initial rates, learn the
//// rates from random initialization
// [[Rcpp::export]]
Rcpp::List nonhomoPois_est_init(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    int H,
    double window,
    double T,
    double dT,
    double gravity,
    arma::cube MuA_start,
    arma::mat tau_init,
    double start,
    bool is_elbo = false
){
  unordered_map<string, std::deque<double>> datamap;
  datamap = transfer_create_empty();
  
  // initialization
  arma::rowvec Pi(K);
  Pi.fill(1.0 / K);
  arma::mat S(m,K);
  arma::cube MuA(K,K,H);
  arma::mat tau(m,K);
  
  // tau.fill(1.0/K);
  tau = tau_init;
  tau = update_init_Tau(tau, m, K);
  for(int k=0; k<K; ++k){
    arma::colvec temp = tau.col(k);
    Pi(k) = mean(temp);
  }
  MuA = MuA_start;
  
  int nall = alltimes.n_rows;
  int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
  // int N = floor(T / dT);
  int N = int((T - start)/dT);
  
  
  arma::vec elbo_vec(N);
  double elbo = 0;
  arma::mat prevdata;
  
  double Tn, t_current, t_start, eta;
  arma::rowvec event; 
  arma::mat truncdata;
  Rcpp::List paralist;
  
  double R = dT;
  
  for (int n = 0; n < N; n++ ){
    Tn = (n + 1.0) * dT;
    event = alltimes.row(start_pos);
    t_current = event(2);
    while (t_current <= Tn ) {
      if (curr_pos >= nall - 1) {
        break;
      } else {
        curr_pos += 1;
        event = alltimes.row(curr_pos);
        t_current = event(2);
      }
    }
    end_pos = curr_pos;
    
    if (end_pos <= start_pos)
      continue;
    
    truncdata = alltimes.rows(start_pos, end_pos - 1);
    
    // datamap = transfer_eff(datamap, truncdata, R);
    transfer_dynamic(datamap, truncdata, R, Tn);
    
    t_start = Tn - dT;
    ln_curr = end_pos;
    n_t = ln_curr - ln_prev;
    eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
    // paralist = update_nonhomo_pois(tau, MuA, Pi, S, datamap, t_start, Tn, m, K, A, window, eta, gravity);
    paralist = update_nonhomo_pois_revised(tau, MuA, Pi, S, datamap, t_start, Tn, m, K, A, window, eta, gravity);
    arma::mat tau_new = paralist["tau"], S_new = paralist["S"];
    arma::cube MuA_new = paralist["MuA"];
    arma::rowvec Pi_new = paralist["Pi"];
    tau = tau_new; 
    // to prevent degenerate groups
    for(int i =0; i < m; ++i){
      tau.row(i) = correct_tau(tau.row(i));
    }
    MuA = MuA_new, S = S_new, Pi = Pi_new;
    start_pos = curr_pos;
    ln_prev = ln_curr;
    // Rprintf("iter: %d; number: %d \n", n, n_t); 
    //MuA.print();
    //S.print();
    
    if (is_elbo){
      prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
      // elbo = get_elbo_nonhomoHak(prevdata, 0, T, tau, MuA, B, Pi, A, lam, m, K, H, window);
      elbo = 0;
      elbo_vec(n) = elbo / ln_curr;
    }
    Rprintf("=============\n");
  }
  
  return Rcpp::List::create(
    Rcpp::Named("MuA") = MuA,
    Rcpp::Named("Pi") = Pi,
    Rcpp::Named("tau") = tau,
    Rcpp::Named("elbo") = elbo_vec);
}
