#include "onlineblock.h"



// [[Rcpp::export]]
arma::vec sampleHak_pre(
    double tStart,
    double T,
    double mu,
    double b,
    double lam
){
  arma::vec parents = genpois(tStart, T, mu); // only change needed
  arma::vec pts_all = parents;
  int loop = 0;
  int N = parents.n_elem;
  if (N == 0) {
    return pts_all;
  } else {
    // printf("loop: %d, number: %d \n", loop, N);
    arma::vec pts_cur = parents;
    int n_curr = pts_cur.n_elem;
    double tcurr;
    while (n_curr > 0){
      arma::vec pts_new, temp;
      for (int i = 0; i < n_curr; i++){
        tcurr = pts_cur(i);
        temp = sampleChild(tcurr,T,lam,b);
        pts_new = vecadd(pts_new, temp);
      }
      pts_all = vecadd(pts_all, pts_new);
      pts_cur = pts_new;
      loop++;
      n_curr = pts_cur.n_elem;
      // printf("loop: %d, number: %d \n", loop, n_curr);
    }
    return arma::sort(pts_all);
  }
}

// generate for all nodes
// want to return the events along with a count matrix of events

// [[Rcpp::export]]
arma::mat sampleBlockHak_pre(
    double T,
    double startT,
    Rcpp::List A,
    arma::vec Z,
    arma::mat Mu,
    arma::mat B,
    double lam
){
  int count = 0, size = 2;
  arma::mat alltimes(size, 3);
  int m = A.size();
  printf("m: %d", m);
  int z1, z2, n_edge, n_temp;
  double mu, b;
  int i,j,k,p;
  arma::vec temp;
  for (i = 0; i < m; i++){
    //printf("count: %d", count);
    arma::rowvec edge = A[i];
    n_edge = edge.n_elem;
    for (p = 0; p < n_edge; p++) {
      j = edge(p);
      z1 = Z(i), z2 = Z(j);
      mu = Mu(z1, z2), b = B(z1, z2);
      temp = sampleHak_pre(startT,T, mu, b, lam);
      n_temp = temp.n_elem;
      //printf("n_temp: %d", n_temp);
      if(count + n_temp > size){
        size = max(2*size, count + n_temp);
        alltimes.resize(size,3);
      }
      for (k = 0; k < n_temp; k++){
        alltimes(count+k,0) = i;
        alltimes(count+k,1) = j;
        alltimes(count+k,2) = temp(k);
      }
      count += n_temp;
    }
  }
  alltimes = alltimes.head_rows(count); // throw away useless rows
  arma::mat alltimes_sort = alltimes; // sort by the third column
  arma::uvec indices = arma::sort_index(alltimes.col(2));
  for (i = 0; i < count; i++) {
    alltimes_sort(i,0) = alltimes(indices(i),0);
    alltimes_sort(i,1) = alltimes(indices(i),1);
    alltimes_sort(i,2) = alltimes(indices(i),2);
  }
  return alltimes_sort;
  
}



// [[Rcpp::export]]
arma::vec sampleHak_nonhomo_pre(
    double T,
    double startT,
    arma::vec avec,
    double b,
    double window,
    double lam
){
  arma::vec parents = genpois_nonhomo(startT, T, window, avec); // only change parents
  arma::vec pts_all = parents;
  int loop = 0;
  int N = parents.n_elem;
  if (N == 0) {
    return pts_all;
  } else {
    // printf("loop: %d, number: %d \n", loop, N);
    arma::vec pts_cur = parents;
    int n_curr = pts_cur.n_elem;
    double tcurr;
    while (n_curr > 0){
      arma::vec pts_new, temp;
      for (int i = 0; i < n_curr; i++){
        tcurr = pts_cur(i);
        temp = sampleChild(tcurr,T,lam,b);
        pts_new = vecadd(pts_new, temp);
      }
      pts_all = vecadd(pts_all, pts_new);
      pts_cur = pts_new;
      loop++;
      n_curr = pts_cur.n_elem;
      // printf("loop: %d, number: %d \n", loop, n_curr);
    }
    return arma::sort(pts_all);
  }
}


// [[Rcpp::export]]
arma::mat sampleBlockHak_nonhomo_pre(
    double T,
    double startT,
    Rcpp::List A,
    arma::vec Z,
    arma::cube MuA,
    arma::mat B,
    double window,
    double lam
){
  int count = 0, size = 2;
  arma::mat alltimes(size, 3);
  int m = A.size(), H = MuA.n_slices;
  printf("H: %d", H);
  printf("m: %d", m);
  int z1, z2, n_edge, n_temp;
  double b;
  int i,j,k,p,h;
  arma::vec temp, avec(H);
  for (i = 0; i < m; i++){
    //printf("count: %d", count);
    arma::rowvec edge = A[i];
    n_edge = edge.n_elem;
    for (p = 0; p < n_edge; p++) {
      j = edge(p);
      z1 = Z(i), z2 = Z(j);
      b = B(z1, z2);
      for (h = 0; h < H; h++) {
        avec(h) = MuA(z1,z2,h);
      }
      temp = sampleHak_nonhomo_pre(T, startT, avec, b, window, lam);
      n_temp = temp.n_elem;
      //printf("n_temp: %d", n_temp);
      if(count + n_temp > size){
        size = max(2*size, count + n_temp);
        alltimes.resize(size,3);
      }
      for (k = 0; k < n_temp; k++){
        alltimes(count+k,0) = i;
        alltimes(count+k,1) = j;
        alltimes(count+k,2) = temp(k);
      }
      count += n_temp;
    }
  }
  alltimes = alltimes.head_rows(count); // throw away useless rows
  arma::mat alltimes_sort = alltimes; // sort by the third column
  arma::uvec indices = arma::sort_index(alltimes.col(2));
  for (i = 0; i < count; i++) {
    alltimes_sort(i,0) = alltimes(indices(i),0);
    alltimes_sort(i,1) = alltimes(indices(i),1);
    alltimes_sort(i,2) = alltimes(indices(i),2);
  }
  return alltimes_sort;
}

