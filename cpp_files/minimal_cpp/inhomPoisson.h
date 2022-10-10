#include "onlineblock.h"


double ELBO_nonhomoHak(
    unordered_map<string, arma::vec> datamap,
    double t_start,
    double t_end,
    arma::mat tau,
    arma::cube MuA,
    arma::mat B,
    arma::rowvec Pi,
    Rcpp::List A,
    double lam,
    int m,
    int K,
    int H,
    double window
){
  double elbo = 0.0;
  
  int l, k, n_edge;
  arma::vec timevec;
  string key;
  int ln,n,h;
  double intensity;
  double period = H * window;
  unordered_map<string, arma::vec>:: iterator itr; 
  
  arma::mat part1(K,K), part2(K,K), Lambda(K,K);
  
  for (itr = datamap.begin(); itr != datamap.end(); itr++) {
    key = itr->first;
    timevec = itr->second;
    arma::vec index = split(key);
    int i = (int) index(0), j = (int) index(1);
    ln = timevec.n_elem;
    
    part1.fill(0.0), part2.fill(0.0);
    for (n = ln - 1; n >= 0; n--){
      double t_current = timevec(n);
      h = floor((t_current - floor(t_current/period) * period)/window);
      if (t_current > t_start) {
        intensity = 0.0;
        Lambda.fill(eps); // store temporary intensity values
        for (int n1 = 0; n1 < n; n1++) {
          double t1 = timevec(n1);
          intensity += trigger(t1, t_current, lam);
        }
        for (k = 0; k < K; k++){
          for (l = 0; l < K; l++){
            Lambda(k,l) += MuA(k,l,h) + B(k,l) * intensity;
          }
        }
        part1 = part1 + arma::log(Lambda);
        part2 = part2 + B * integral(t_current, t_end, lam);
      } else {
        part2 = part2 + B * integral2(t_current, t_start, t_end, lam);
      }
    }
    
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        elbo += tau(i,k) * tau(j,l) * (part1(k,l) - part2(k,l));
      }
    }
  }
  
  // mu * T，
  
  arma::vec tvec(H);
  tvec.fill(0.0);
  int h1 = floor(t_start/window);
  int h2 = floor(t_end/window);
  
  for (int w = h1; w <= h2; w++){
    h = w % H;
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
  
  for (int i = 0; i < m; i++) {
    arma::rowvec edge = A[i];
    n_edge = edge.n_elem;
    for (int p = 0; p < n_edge; p++) {
      int j = (int) edge(p);
      for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
          for (h = 0; h < H; h++)
            elbo -= tau(i,k) * tau(j,l) * MuA(k,l,h) * tvec(h);
        }
      }
    }
    
  }
  
  // tau
  for (int i = 0; i < m; i++) {
    for (k = 0; k < K; k++) {
      elbo += tau(i,k) * (log(tau(i,k) + eps) - log(Pi(k) + eps));
    }
  }
  
  return elbo;
}


// [[Rcpp::export]]
double get_elbo_nonhomoHak(
    arma::mat alltimes,
    double t_start,
    double t_end,
    arma::mat tau,
    arma::cube MuA,
    arma::mat B,
    arma::rowvec Pi,
    Rcpp::List A,
    double lam,
    int m,
    int K,
    int H,
    double window
){
  unordered_map<string, arma::vec> datamap = transfer(alltimes);
  double elbo = ELBO_nonhomoHak(datamap, t_start, t_end, tau, MuA, B, Pi, A, lam, m, K, H, window);
  return elbo;
}



// poisson case

Rcpp::List update_nonhomo_pois(
    arma::mat tau,
    arma::cube MuA,
    arma::rowvec Pi,
    arma::mat S,
    unordered_map<string, std::deque<double>> datamap,
    double t_start,
    double t_end,
    int m,
    int K,
    Rcpp::List A,
    double window,
    double eta,
    double gravity
){
  int H = MuA.n_slices;
  double period = H * window;
  arma::cube P1_mu(K,K,H), P2_mu(K,K,H), P1_mu_tp(K,K,H), P2_mu_tp(K,K,H);
  arma::mat S_tp(m,K);
  arma::mat P_S_tp(m,K), Lambda(K,K);
  P1_mu.fill(0.0), P2_mu.fill(0.0); 
  S_tp.fill(0.0), Lambda.fill(0.0);
  
  int l, k, h;
  arma::vec tvec(H);
  tvec.fill(0.0);
  int h1 = floor(t_start/window);
  int h2 = floor(t_end/window);
  
  //Rprintf("h1: %d,  h2: %d \n", h1, h2);
  for (int w = h1; w <= h2; w++){
    h = w % H;
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
  
  /*
   for (int i = 0; i < m; i++) {
   arma::rowvec edge = A[i];
   n_edge = edge.n_elem;
   for (int p = 0; p < n_edge; p++){
   int j = (int) edge(p);
   for (k = 0; k < K; k++) {
   for (l = 0; l < K; l++) {
   for (h = 0; h < H; h++)
   P2_mu(k,l,h) = P2_mu(k,l,h) + tau(i,k) * tau(j,l) * tvec(h);
   }
   }
   }
   }
   */
  
  unordered_map<string, std::deque<double>>:: iterator itr; 
  std::deque<double> timeque;
  arma::vec timevec;
  string key;
  int ln,n;
  
  
  arma::cube MuA_count(K,K,H);
  MuA_count.fill(0.0);
  int n_events = 0;
  
  for (itr = datamap.begin(); itr != datamap.end(); itr++) 
  { 
    // type itr->first stores the key part  and 
    // itr->second stroes the value part 
    key = itr->first;
    timeque = itr->second;
    timevec = convert_deque(timeque);
    
    arma::vec index = split(key);
    int i = (int) index(0), j = (int) index(1);
    
    if (i == j)
      continue;
    
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        for (h = 0; h < H; h++)
          P2_mu(k,l,h) = P2_mu(k,l,h) + tau(i,k) * tau(j,l) * tvec(h);
      }
    }
    
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        for (h = 0; h < H; h++) {
          S_tp(i,k) = S_tp(i,k) - tau(j,l) * MuA(k,l,h) * tvec(h); // change here
        }
      }
    }
    
    P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0); 
    P_S_tp.fill(0.0);
    
    ln = timevec.n_elem;
    
    for (n = ln - 1; n >= 0; n--){
      double t_current = timevec(n);
      h = floor((t_current - floor(t_current/period) * period)/window); // which window
      
      for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++)
          MuA_count(k,l,h) += 1;
      }
      n_events += 1;
      
      if (t_current > t_start) {
        Lambda.fill(eps); // store temporary intensity values
        
        for (k = 0; k < K; k++){
          for (l = 0; l < K; l++){
            Lambda(k,l) += MuA(k,l,h); // changed 
          }
        }
        
        for (k = 0; k < K; k++) {
          for (l = 0; l < K; l++){
            P1_mu_tp(k,l,h) += 1.0 / Lambda(k,l);
            P_S_tp(k,l) += log(Lambda(k,l));
          }
        }
      } else {
        break;
      }
    }
    
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        for (h = 0; h < H; h++) {
          P1_mu(k,l,h) += tau(i,k) * tau(j,l) * P1_mu_tp(k,l,h);
        }
      }
    }
    
    // update S
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        S_tp(i,k) += tau(j,l) * (P_S_tp(k,l)); // change here
      }
    }
  } 
  
  /*
   // update S, second part
   for (int i = 0; i < m; i++) {
   arma::rowvec edge = A[i];
   n_edge = edge.n_elem;
   for (k = 0; k < K; k++) {
   for (int p = 0; p < n_edge; p++) {
   int j = (int) edge(p);
   for (l = 0; l < K; l++) {
   for (h = 0; h < H; h++) {
   S_tp(i,k) = S_tp(i,k) - tau(j,l) * MuA(k,l,h) * tvec(h); // change here
   }
   }
   }
   }
   }
   */
  
  // update parameters
  S = S + S_tp;
  arma::cube grad_mu = P1_mu - P2_mu;
  //grad_mu.print();
  arma::cube MuA_new = MuA + eta * grad_mu;
  //Rprintf("Mu new is: \n");
  //Mu_new.print();
  
  // handle negative values and large gradient
  for (k = 0; k < K; k++) {
    for (l = 0; l < K; l++) {
      for (h = 0; h < H; h++) {
        if (MuA_new(k,l,h) <= 0.0)
          MuA_new(k,l,h) = MuA(k,l,h) / 2.0;
        else if (MuA_new(k,l,h) > 2 * MuA(k,l,h))
          MuA_new(k,l,h) = MuA(k,l,h) * 2.0;
        
        if (MuA_new(k,l,h) <= gravity * MuA_count(k,l,h) / n_events)
          MuA_new(k,l,h) = 0.0;
      }
    }
  }
  
  
  arma::mat tau_new(m,K);
  tau_new.fill(0.0);
  
  arma::rowvec s;
  for (int i = 0; i < m; i++) {
    arma::rowvec s = arma::log(Pi) + S.row(i);
    s = s - max(s);
    s = exp(s)/sum(exp(s));
    tau_new.row(i) = correct_tau(s);
  }
  
  for (k = 0; k < K; k++) {
    Pi(k) = sum(tau_new.col(k)) / (m + 0.0);
  }	
  
  return Rcpp::List::create(Rcpp::Named("tau") = tau_new,
                            Rcpp::Named("MuA") = MuA_new,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("S") = S);
}


Rcpp::List update_nonhomo_pois_revised(
    arma::mat tau,
    arma::cube MuA,
    arma::rowvec Pi,
    arma::mat S,
    unordered_map<string, std::deque<double>> datamap,
    double t_start,
    double t_end,
    int m,
    int K,
    Rcpp::List A,
    double window,
    double eta,
    double gravity
){
  int H = MuA.n_slices;
  double period = H * window;
  arma::cube P1_mu(K,K,H), P2_mu(K,K,H), P1_mu_tp(K,K,H), P2_mu_tp(K,K,H);
  arma::mat S_tp(m,K);
  arma::mat P_S_tp(m,K), Lambda(K,K);
  P1_mu.fill(0.0), P2_mu.fill(0.0); 
  S_tp.fill(0.0), Lambda.fill(0.0);
  
  int l, k, h;
  arma::vec tvec(H);
  tvec.fill(0.0);
  int h1 = floor(t_start/window);
  int h2 = floor(t_end/window);
  
  //Rprintf("h1: %d,  h2: %d \n", h1, h2);
  for (int w = h1; w <= h2; w++){
    h = w % H;
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
  
  int n_edge;
  for (int i = 0; i < m; i++) {
    arma::rowvec edge = A[i];
    n_edge = edge.n_elem;
    for (int p = 0; p < n_edge; p++){
      int j = (int) edge(p);
      if(i != j){
        for (k = 0; k < K; k++) {
          for (l = 0; l < K; l++) {
            for (h = 0; h < H; h++)
              P2_mu(k,l,h) = P2_mu(k,l,h) + tau(i,k) * tau(j,l) * tvec(h);
          }
        }
        
        for (k = 0; k < K; k++) {
          for (l = 0; l < K; l++) {
            for (h = 0; h < H; h++) {
              S_tp(i,k) = S_tp(i,k) - tau(j,l) * MuA(k,l,h) * tvec(h); // change here
            }
          }
        }
      }
    }
    
  }
  
  
  unordered_map<string, std::deque<double>>:: iterator itr; 
  std::deque<double> timeque;
  arma::vec timevec;
  string key;
  int ln,n;
  
  
  arma::cube MuA_count(K,K,H);
  MuA_count.fill(0.0);
  int n_events = 0;
  
  for (itr = datamap.begin(); itr != datamap.end(); itr++) 
  { 
    // type itr->first stores the key part  and 
    // itr->second stroes the value part 
    key = itr->first;
    timeque = itr->second;
    timevec = convert_deque(timeque);
    
    arma::vec index = split(key);
    int i = (int) index(0), j = (int) index(1);
    
    if (i == j)
      continue;
    
    /*
     for (k = 0; k < K; k++) {
     for (l = 0; l < K; l++) {
     for (h = 0; h < H; h++)
     P2_mu(k,l,h) = P2_mu(k,l,h) + tau(i,k) * tau(j,l) * tvec(h);
     }
     }
     
     for (k = 0; k < K; k++) {
     for (l = 0; l < K; l++) {
     for (h = 0; h < H; h++) {
     S_tp(i,k) = S_tp(i,k) - tau(j,l) * MuA(k,l,h) * tvec(h); // change here
     }
     }
     }
     */
    
    P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0); 
    P_S_tp.fill(0.0);
    
    ln = timevec.n_elem;
    
    for (n = ln - 1; n >= 0; n--){
      double t_current = timevec(n);
      h = floor((t_current - floor(t_current/period) * period)/window); // which window
      
      for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++)
          MuA_count(k,l,h) += 1;
      }
      n_events += 1;
      
      if (t_current > t_start) {
        Lambda.fill(eps); // store temporary intensity values
        
        for (k = 0; k < K; k++){
          for (l = 0; l < K; l++){
            Lambda(k,l) += MuA(k,l,h); // changed 
          }
        }
        
        for (k = 0; k < K; k++) {
          for (l = 0; l < K; l++){
            P1_mu_tp(k,l,h) += 1.0 / Lambda(k,l);
            P_S_tp(k,l) += log(Lambda(k,l));
          }
        }
      } else {
        break;
      }
    }
    
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        for (h = 0; h < H; h++) {
          P1_mu(k,l,h) += tau(i,k) * tau(j,l) * P1_mu_tp(k,l,h);
        }
      }
    }
    
    // update S
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        S_tp(i,k) += tau(j,l) * (P_S_tp(k,l)); // change here
      }
    }
  } 
  
  /*
   // update S, second part
   for (int i = 0; i < m; i++) {
   arma::rowvec edge = A[i];
   n_edge = edge.n_elem;
   for (k = 0; k < K; k++) {
   for (int p = 0; p < n_edge; p++) {
   int j = (int) edge(p);
   for (l = 0; l < K; l++) {
   for (h = 0; h < H; h++) {
   S_tp(i,k) = S_tp(i,k) - tau(j,l) * MuA(k,l,h) * tvec(h); // change here
   }
   }
   }
   }
   }
   */
  
  // update parameters
  S = S + S_tp;
  arma::cube grad_mu = P1_mu - P2_mu;
  //grad_mu.print();
  arma::cube MuA_new = MuA + eta * grad_mu;
  //Rprintf("Mu new is: \n");
  //Mu_new.print();
  
  // handle negative values and large gradient
  for (k = 0; k < K; k++) {
    for (l = 0; l < K; l++) {
      for (h = 0; h < H; h++) {
        if (MuA_new(k,l,h) <= 0.0)
          MuA_new(k,l,h) = MuA(k,l,h) / 2.0;
        else if (MuA_new(k,l,h) > 2 * MuA(k,l,h))
          MuA_new(k,l,h) = MuA(k,l,h) * 2.0;
        
        if (MuA_new(k,l,h) <= gravity * MuA_count(k,l,h) / n_events)
          MuA_new(k,l,h) = 0.0;
      }
    }
  }
  
  
  arma::mat tau_new(m,K);
  tau_new.fill(0.0);
  
  arma::rowvec s;
  for (int i = 0; i < m; i++) {
    arma::rowvec s = arma::log(Pi) + S.row(i);
    s = s - max(s);
    s = exp(s)/sum(exp(s));
    tau_new.row(i) = correct_tau(s);
  }
  
  for (k = 0; k < K; k++) {
    Pi(k) = sum(tau_new.col(k)) / (m + 0.0);
  }	
  
  return Rcpp::List::create(Rcpp::Named("tau") = tau_new,
                            Rcpp::Named("MuA") = MuA_new,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("S") = S);
}



// [[Rcpp::export]]
Rcpp::List nonhomoPois_estimator(
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
    arma::mat tau_start,
    bool is_elbo = false
){
  unordered_map<string, std::deque<double>> datamap;
  datamap = transfer_create_empty();
  
  // initialization
  arma::rowvec Pi(K);
  Pi.fill(1.0 / K);
  // Pi = {0.2, 0.8};
  arma::mat S(m,K);
  arma::cube MuA(K,K,H);
  arma::mat tau(m,K);
  S.fill(1.0/K);
  for (int k = 0; k < K; k++) {
    for (int l=0; l < K; l++) {
      for(int h = 0; h < H; h++) {
        MuA(k,l,h) = myrunif();
      }
    }
  }
  tau.fill(1.0/K);
  //tau = tau_start;
  
  int nall = alltimes.n_rows;
  int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
  int N = floor(T / dT);
  
  
  arma::vec elbo_vec(N);
  double elbo = 0;
  arma::mat prevdata;
  
  double Tn, t_current, t_start, eta;
  arma::rowvec event; 
  arma::mat truncdata;
  Rcpp::List paralist;
  
  double R = dT;
  
  for (int n = 0; n < N; n++ ){
    cout<<n<<endl;
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
    MuA = MuA_new, S = S_new, Pi = Pi_new;
    start_pos = curr_pos;
    ln_prev = ln_curr;
    // Rprintf("iter: %d; number: %d \n", n, n_t); 
    // MuA.print();
    // Pi.print();
    
    if (is_elbo){
      // cumulative elbo for all previous data
      prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
      // need to declare B and lam as 0 in here
      arma::mat B(K, K);
      B.fill(0);
      double lam = 0;
      elbo = get_elbo_nonhomoHak(prevdata, 0, T, tau, MuA, B, Pi, A, lam, m, K, H, window);
      // elbo = 0;
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

//// modify inhomogeneous Poisson to take initial values also
//// only pass in the community assignments, learn the
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
  MuA.fill(0.5), S.fill(0.0);
  for (int k = 0; k < K; k++) {
    for (int l=0; l < K; l++) {
      for(int h = 0; h < H; h++) {
        MuA(k,l,h) = myrunif();
      }            
    }
  }
  // tau.fill(1.0/K);
  tau = tau_init;
  
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
    MuA = MuA_new, S = S_new, Pi = Pi_new;
    start_pos = curr_pos;
    ln_prev = ln_curr;
    Rprintf("iter: %d; number: %d \n", n, n_t); 
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


// [[Rcpp::export]]
Rcpp::List batch_nonhomoPois_estimator(
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
    arma::mat tau_start,
    int itermax,
    double stop_eps
){
  // initialization
  arma::rowvec Pi(K);
  Pi.fill(1.0 / K);
  arma::mat S(m,K);
  arma::cube MuA(K,K,H);
  arma::mat tau(m,K);
  MuA.fill(0.5), S.fill(0.0);
  
  for (int k = 0; k < K; k++) {
    for (int l=0; l < K; l++) {
      for (int h=0; h <H ; h++)
        MuA(k,l,h) = myrunif();
    }
  }
  //B = B_start, Mu = Mu_start;
  // for (int i = 0; i < m; i++) {
  // 	arma::rowvec tt(K);
  // 	for (int k = 0; k < K; k++) {
  // 		tt(k) = myrunif();
  // 	}
  // 	tt = tt / sum(tt);
  // 	tau.row(i) = tt;
  // }
  tau.fill(1.0/K);
  //tau = tau_start;
  
  int nall = alltimes.n_rows;
  int ncol = alltimes.n_cols;
  
  Rcpp::List paralist;
  
  unordered_map<string, std::deque<double>> datamap;
  datamap = transfer_create(A,m);
  double R = T;
  // datamap = transfer_eff(datamap, alltimes, R);
  transfer_eff(datamap, alltimes, R);
  
  double t_start = 0.0, Tn;
  Tn = alltimes(nall - 1, ncol - 1);
  
  double gap = 2147483647;
  double eta = 1.0/nall * (K * K);
  
  for (int iter = 0; iter < itermax; iter++) {
    eta = 1.0/nall * (K * K) / sqrt(iter / H + 1.0);
    S.fill(0.0);        
    paralist = update_nonhomo_pois(tau, MuA, Pi, S, datamap, t_start, Tn, m, K, A, window, eta, gravity);
    // paralist = update_nonhomo_sparse(tau, MuA, B, Pi, S, datamap, t_start, Tn, m, K, A, window, lam, eta, gravity);
    arma::mat tau_new = paralist["tau"], S_new = paralist["S"];
    arma::cube MuA_new = paralist["MuA"];
    arma::rowvec Pi_new = paralist["Pi"];
    
    gap = abs(MuA - MuA_new).max();
    tau = tau_new; 
    MuA = MuA_new, S = S_new, Pi = Pi_new;
    Rprintf("iter: %d \n", iter); 
    // B.print();
    // Mu.print();
    Pi.print();
    Rprintf("gap: %2.3f", gap);
    Rprintf("=============\n");
    if (gap < stop_eps){
      break;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("MuA") = MuA,
    Rcpp::Named("Pi") = Pi,
    Rcpp::Named("tau") = tau);
}