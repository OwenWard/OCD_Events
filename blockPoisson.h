#include "onlineblock.h"




// [[Rcpp::export]]
arma::rowvec initPi(int K){
  arma::rowvec Pi(K);
  Pi.fill(1.0);
  return Pi/K;
}


// [[Rcpp::export]]
arma::mat initS(int m, int K){
  arma::mat S(m,K);
  return S.randu(m,K);
}

// [[Rcpp::export]]
arma::mat initB(int K){
  arma::mat B(K,K);
  return B.randu(K,K);
}




// [[Rcpp::export]]
arma::mat updateTau(arma::mat S, arma::rowvec Pi, int m, int K){
  arma::mat tau;
  tau.zeros(m,K);
  for(int i =0; i < m; ++i){
    arma::rowvec s = arma::log(Pi) + S.row(i);
    s = s - max(s);
    s = exp(s)/sum(exp(s));
    tau.row(i) = s;
  }
  return tau;
}

// [[Rcpp::export]]
arma::rowvec updatePi(arma::mat tau, int K){
  arma::rowvec Pi(K);
  for(int k=0; k<K; ++k){
    arma::colvec temp = tau.col(k);
    Pi(k) = mean(temp);
  }
  return Pi;
}

// [[Rcpp::export]]
arma::mat updateS(arma::mat data, arma::mat tau, arma::mat B,
                  Rcpp::List A, //arma::mat A,
                  arma::mat S, int K, int m, double dT){
  int nrow = data.n_rows;
  for(int n = 0; n < nrow; ++n){
    arma::rowvec event = data.row(n);
    int i = event(0);
    int j = event(1);
    for(int k = 0; k < K; ++k){
      for(int l = 0; l < K; ++l){
        S(i,k) = S(i,k) + tau(j,l)*log(B(k,l)); 
      }
    }
  }
  for(int i=0; i < m; ++i){
    //arma::rowvec edge = A.row(i);
    arma::rowvec edge = A[i];
    int n = edge.n_elem;
    for(int k = 0; k < K; ++k){
      for(int j = 0; j < n; ++j){
        for(int l =0; l<K; ++l){
          int j_loc = edge(j);
          S(i,k) = S(i,k) - tau(j_loc,l)*B(k,l)*dT;
        }
      }
    }
  }
  return S;
}

// [[Rcpp::export]]
arma::mat updateB(arma::mat data, arma::mat tau, arma::mat B, int K,
                  Rcpp::List A, //arma::mat A,
                  int m , double dT, double eta){
  arma::mat X1,X2, gradB, B_new;
  X1.zeros(K,K);
  X2.zeros(K,K);
  gradB.zeros(K,K);
  int nrow = data.n_rows;
  for(int n = 0; n < nrow; n++){
    arma::rowvec event = data.row(n);
    int i = event(0);
    int j = event(1);
    for(int k =0; k<K; ++k){
      for(int l =0; l<K; ++l){
        X1(k,l) += tau(i,k)*tau(j,l);
      }
    }
  }
  //auto start = high_resolution_clock::now();
  for(int i =0; i<m; ++i){
    //arma::rowvec edge = A.row(i);
    arma::rowvec edge = A[i];
    int n = edge.n_elem;
    for(int j=0; j< n; ++j){
      int j_loc = edge(j);
      //cout<<i<<" is i index "<<j_loc<<" is j index"<<endl;
      if( i != j_loc){
        for(int k=0; k<K; ++k){
          for(int l =0; l<K; ++l){
            X2(k,l) += tau(i,k)*tau(j_loc,l)*dT;
          }
        }
      }
    }
  }
  // auto stop = high_resolution_clock::now();
  // auto duration = duration_cast<milliseconds>(stop-start);
  // cout<<duration.count()<<endl;
  gradB = X1/B - X2;
  //cout << gradB << endl;
  B_new = B + eta*gradB; //previously had division of n_row here
  // create matrix of 0.001 then take max element wise
  // prevent large updates
  for (int k = 0; k < K; k++) {
    for (int l = 0; l < K; l++) {
      if (B_new(k,l) <= 0.0)
        B_new(k,l) = B(k,l) / 2.0;
      if (B_new(k,l) > 2 * B(k,l))
        B_new(k,l) = B(k,l) * 2.0;
    }
  }
  //arma::mat eps(K,K);
  //eps.fill(0.001);
  // B_new = arma::max(eps,B_new);
  return B_new;
}


// then the full thing in one function

// [[Rcpp::export]]
Rcpp::List update_Poisson(
    arma::mat data, arma::mat tau, 
    arma::mat B, arma::rowvec Pi,
    arma::mat S, Rcpp::List A, //arma::mat A,
    int m, int K, double eta, double dT){
  
  S = updateS(data,tau,B,A,S,K,m,dT);
  tau = updateTau(S,Pi,m,K); 
  B = updateB(data,tau,B,K,A,m,dT,eta);
  Pi = updatePi(tau,K);
  
  
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            Named("B")=B,
                            Named("Pi")=Pi);
}


// [[Rcpp::export]]
double computeELBO(
    arma::mat data, arma::mat tau,
    arma::mat B, arma::rowvec Pi, Rcpp::List A,
    int m, int K, double dT){
  double epsilon = 0.0001;
  double elbo = 0;
  int nrow = data.n_rows;
  for(int n =0; n<nrow; n++){
    arma::rowvec event = data.row(n);
    int i = event(0);
    int j = event(1);
    for(int k = 0; k <K; ++k){
      for(int l =0; l<K; l++){
        elbo += tau(i,k)*tau(j,l)*log(B(k,l));
      }
    }
  }
  for(int i=0; i<m; ++i){
    arma::rowvec edge = A[i];
    int n = edge.n_elem;
    for(int k=0; k< K; ++k){
      for(int j=0; j<n; ++j){
        int j_loc = edge(j);
        if(i != j_loc){
          for(int l =0; l<K; ++l){
            elbo += -tau(i,k)*tau(j_loc,l)*B(k,l)*dT;
          }
        }
      }
      elbo += tau(i,k)*(log(Pi(k))-log(tau(i,k)+epsilon));
    }
  }
  return elbo;
}

// [[Rcpp::export]]
double computeLL(arma::mat data, arma::mat tau,
                 arma::mat B, arma::rowvec Pi, Rcpp::List A,
                 int m, int K, double currT){
  // compute the log likelihood given current estimates
  double LL = 0;
  arma:: vec curr_group;
  curr_group.zeros(m);
  for(int i = 0; i < m; ++i){
    arma::rowvec curr_tau = tau.row(i);
    curr_group(i) = curr_tau.index_max(); // need this later anyway
    LL += log(Pi(curr_group(i)));
  }
  // compute the next term in the same way as done for the elbo.
  int nrow = data.n_rows;
  for(int n =0; n<nrow; n++){
    arma::rowvec event = data.row(n);
    int i = event(0);
    int j = event(1);
    int k = curr_group(i);
    int l = curr_group(j);
    LL += log(B(k,l));
  }
  for(int i=0; i<m; ++i){
    arma::rowvec edge = A[i];
    int n = edge.n_elem;
    for(int j=0; j<n; ++j){
      int j_loc = edge(j);
      if(i != j_loc){
        int k = curr_group(i);
        int l = curr_group(j_loc);
        LL -= B(k,l)*currT;
      
      }
    }
  }
  return LL;
}


// [[Rcpp::export]]
Rcpp::List estimate_Poisson(
    arma::mat full_data, 
    Rcpp::List A, //arma::mat A,
    int m,
    int K,
    double T,
    double dT,
    arma::mat B,
    arma::mat tau, 
    arma::rowvec Pi,
    arma::mat S,
    int inter_T,
    bool is_elbo = false
  ){
  // iterate this over the time windows...
  int N = int(T/dT);
  int slices = int(N/inter_T);
  double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  int ind = 0;
  int nall = full_data.n_rows;
  arma::cube inter_tau(m,K,slices+1);
  arma::cube inter_B(K,K,N);
  arma::vec curr_elbo, ave_elbo, ave_ll, curr_ll;
  curr_elbo.zeros(N);
  curr_ll.zeros(N);
  ave_ll.zeros(N);
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
    sub_data = full_data.rows(start_pos,end_pos);
    cum_events += sub_data.n_rows;
    elbo_dat = full_data.rows(0,end_pos); 
    //cout<<size(sub_data)<<endl;
    start_pos = curr_pos;
    eta = 1/pow(1+n, .5)/sub_data.n_rows*(K*K);
    S = updateS(sub_data,tau,B,A,S,K,m,dT);
    //cout<<"S works"<<endl;
    tau = updateTau(S,Pi,m,K); 
    //cout<<"update tau"<<endl;
    B = updateB(sub_data,tau,B,K,A,m,dT,eta);
    inter_B.slice(n) = B;
    //cout<<"update B"<<endl;
    Pi = updatePi(tau,K);
    if (is_elbo) {
      curr_elbo(n) = computeELBO(elbo_dat,tau,B,Pi,A,m,K,dT);
      ave_elbo(n) = curr_elbo(n)/cum_events;
      curr_ll(n) = computeLL(elbo_dat,tau,B,Pi,A,m,K,t_curr);
      ave_ll(n) = curr_ll(n)/cum_events;
    }
    
    //cout<<B<<endl;
    //printf("iter: %d; \n", n); 
    //B.print();
    //Pi.print();
    //S.print();
    //printf("=============\n");
    if(n % inter_T == 0 ){
      inter_tau.slice(ind) = tau;
      ind = ind + 1;
      printf("iter: %d; \n", n);
      printf("=============\n");
    }
    
  }
  
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            Named("early_tau")= inter_tau,
                            Named("inter_B") = inter_B,
                            Named("B")=B,
                            Named("Pi")=Pi,
                            Named("AveELBO")=ave_elbo,
                            Named("logL") = ave_ll);
}




// batch optimization
// [[Rcpp::export]]
Rcpp::List batch_estimator_hom_Poisson(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    arma::mat B_start,
    arma::mat tau_start,
    int itermax,
    double stop_eps
){
  // initialization
  arma::rowvec Pi(K);
  Pi.fill(1.0 / K);
  arma::mat B(K,K), S(m,K);
  arma::mat tau(m,K);
  for (int k = 0; k < K; k++) {
    for (int l=0; l < K; l++) {
      B(k,l) = myrunif();
    }
  }
  //B.fill(0.5), Mu.fill(0.5); 
  
  //B = B_start, Mu = Mu_start;
  for (int i = 0; i < m; i++) {
    arma::rowvec tt(K);
    for (int k = 0; k < K; k++) {
      tt(k) = myrunif();
    }
    tt = tt / sum(tt);
    tau.row(i) = tt;
  }
  
  int nall = alltimes.n_rows;
  double gap = 2147483647;
  double eta = 1.0/nall * (K * K);
  
  arma::vec curr_elbo;
  curr_elbo.zeros(itermax);
  double elbo_gap = 2147483647;
  
  
  for (int iter = 0; iter < itermax; iter++) {
    eta = 1.0/ (K * K) /(iter + 1.0);
    //eta = 0.001;
    //printf("eta: %f",eta);
    
    // then do all the updates in here....
    // what happens to dT?
    S.fill(0.0);
    arma::mat S_new = updateS(alltimes,tau,B,A,S,K,m,T);
    //cout<<"S works"<<endl;
    arma::mat tau_new = updateTau(S_new,Pi,m,K); 
    //cout<<"update tau"<<endl;
    arma::mat B_new = updateB(alltimes,tau_new,B,K,A,m,T,eta);
    //cout<<"update B"<<endl;
    arma::rowvec Pi_new = updatePi(tau_new,K);
    curr_elbo(iter) = computeELBO(alltimes,tau_new,B_new,Pi_new,A,m,K,T);
    //curr_ll(n) = computeLL(alltimes,tau,B,Pi,A,m,K,t_curr);
    // need to check this
    // then compare gap...
    gap = abs(B - B_new).max();
    if(iter > 0){
      elbo_gap = abs(curr_elbo(iter) - curr_elbo(iter-1));
    }
    
    // then convert them
    tau = tau_new, B = B_new, Pi = Pi_new, S = S_new;
    printf("gap: %2.3f", gap);
    printf("=============\n");
    B.print();
    if (gap < stop_eps){
      break;
    }
  }
  
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            //Named("early_tau")= inter_tau,
                            //Named("inter_B") = inter_B,
                            Named("B")=B,
                            Named("Pi")=Pi,
                            Named("ELBO")=curr_elbo);
}


