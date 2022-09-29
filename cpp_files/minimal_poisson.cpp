// full online estimation procedure
// minimal implementation of online poisson
// procedure, to quickly compare and implement
// [[Rcpp::export]]
Rcpp::List estimate_Poisson_minimal(
  arma::mat full_data, 
  Rcpp::List A, //arma::mat A,
  int m,
  int K,
  double T,
  double dT,
  arma::mat B,
  bool is_elbo = false
){
  // iterate this over the time windows...
  int N = int(T/dT);
  double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  int nall = full_data.n_rows;
  arma::rowvec Pi(K);
  Pi.fill(1.0 / K);
  arma::mat S(m,K);
  arma::mat tau(m,K);
  S.fill(1.0/K);
  tau.fill(1.0/K);
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll, full_elbo;
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
  ave_elbo.zeros(N);
  full_elbo.zeros(N);
  int cum_events = 0;
  for(int n = 0; n < N; ++n){
    cout<<n<<endl;
    double Tn = dT*(n+1);
    // double T0 = dT*n;
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
    ////
      // arma::vec times = full_data.col(2);
      // arma::uvec valid = find(times < Tn && times > T0);
      // arma::mat sub_data = full_data.rows(valid);
      // cout<<curr_times.n_rows<<endl;
      // then just need to get the start of the indexing also, use that
      // to get the whole subset
      //
        end_pos = curr_pos;
        arma::mat sub_data, elbo_dat;
        sub_data = full_data.rows(start_pos, end_pos);
        cum_events += sub_data.n_rows;
        elbo_dat = full_data.rows(0,end_pos);
        start_pos = curr_pos;
        eta = 1/pow(1 + n, .5)/sub_data.n_rows*(K*K);
        // cout<<dT<<endl;
        S = updateS(sub_data, tau, B, A, S, K, m, dT);
        tau = updateTau(S,Pi,m,K); 
        B = updateB(sub_data, tau, B, K, A, m, dT, eta);
        Pi = updatePi(tau,K);
        
  }
  
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            Named("B")=B,
                            Named("Pi")=Pi,
                            Named("wind_elbo") = curr_elbo,
                            Named("AveELBO")=ave_elbo,
                            Named("wind_ll") = curr_ll,
                            Named("logL") = ave_ll,
                            Named("full_ELBO") = full_elbo);
}
