#include "onlineblock.h"


// general function to compute log likelihood for inhomogeneous processes

// [[Rcpp::export]]
double ll_nonhomo(
    unordered_map<string, arma::vec> datamap,
    double t_start,
    double t_end,
    arma::vec z_est, // computed in R, so not zero indexed
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
  double ll = 0.0;
  
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
    
    k = z_est(i) - 1; // for zero index
    l = z_est(j) - 1;
    ll += (part1(k,l) - part2(k,l));

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
      for (h = 0; h < H; h++){
        k = z_est(i) - 1; // for zero index
        l = z_est(j) - 1;
        ll -= MuA(k,l,h) * tvec(h);
      }
    }
  }
  
  // tau
  for (int i = 0; i < m; i++) {
    k = z_est(i) - 1;
    ll -= log(Pi(k) + eps));
  }
  
  return ll;
}
