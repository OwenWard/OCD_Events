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
  gradB = X1/B - X2;
  //cout << gradB << endl;
  B_new = B + eta*gradB/nrow;
  // create matrix of 0.001 then take max element wise
  // prevent large updates
  // for (int k = 0; k < K; k++) {
  //   for (int l = 0; l < K; l++) {
  //     if (B_new(k,l) <= 0.0) 
  //       B_new(k,l) = B(k,l) / 2.0;
  //     if (B_new(k,l) > 2 * B(k,l))
  //       B_new(k,l) = B(k,l) * 2.0;
  //   }
  // }
  arma::mat eps(K,K);
  eps.fill(0.001);
  B_new = arma::max(eps,B_new);
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
  // compute the negative log likelihood given current estimates
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
    arma::mat full_data, arma::mat tau, 
    arma::mat B, arma::rowvec Pi,
    arma::mat S, Rcpp::List A, //arma::mat A,
    int m, int K, double dT, double T){
  // iterate this over the time windows...
  int N = int(T/dT);
  int slices = int(N/50);
  double eta;
  int start_pos = 0;
  int curr_pos = 0;
  int end_pos = 0;
  int ind = 0;
  int nall = full_data.n_rows;
  arma::cube inter_tau(m,K,slices+1);
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
    eta = 1/sqrt(1+n);
    S = updateS(sub_data,tau,B,A,S,K,m,dT);
    //cout<<"S works"<<endl;
    tau = updateTau(S,Pi,m,K); 
    //cout<<"update tau"<<endl;
    B = updateB(sub_data,tau,B,K,A,m,dT,eta);
    //cout<<"update B"<<endl;
    Pi = updatePi(tau,K);
    curr_elbo(n) = computeELBO(elbo_dat,tau,B,Pi,A,m,K,dT);
    curr_ll(n) = computeLL(elbo_dat,tau,B,Pi,A,m,K,t_curr);
    ave_ll(n) = curr_ll(n)/cum_events;
    ave_elbo(n) = curr_elbo(n)/cum_events;
    //cout<<B<<endl;
    printf("iter: %d; \n", n); 
    //B.print();
    //Pi.print();
    //S.print();
    printf("=============\n");
    if(n % 50 == 0){
      inter_tau.slice(ind) = tau;
      ind = ind + 1;
    }
    
  }
  
  return Rcpp::List::create(Named("S")= S,
                            Named("tau")=tau,
                            Named("early_tau")= inter_tau,
                            Named("B")=B,
                            Named("Pi")=Pi,
                            Named("AveELBO")=ave_elbo,
                            Named("logL") = ave_ll);
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


