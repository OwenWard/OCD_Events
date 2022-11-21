#include "onlineblock.h"


// [[Rcpp::export]]
arma::mat max_tau(arma::mat tau){
  int n = tau.n_rows;
  int m = tau.n_cols;
  arma::mat tau_ind(n, m);
  tau_ind.fill(0.0);
  for(int i = 0; i < n; ++i){
    arma::rowvec curr = tau.row(i);
    int ind = curr.index_max();
    tau_ind(i,ind) = 1;
  }
  return tau_ind;
}


// [[Rcpp::export]]
arma::cube to_cube(arma::vec values, int x, int y, int z) {
  arma::vec v(values);
  arma::cube res((const double*)v.begin(), x, y, z);
  return res;
}


// [[Rcpp::export]]
arma::rowvec to_vec(arma::cube Q){
  // arma::rowvec A = Q(arma::span(0), arma::span(1), arma::span::all);
  arma::vec A = vectorise(Q);
  return A.as_row();
}

// [[Rcpp::export]]
double myrexp(
	double lam
	){
	// random device class instance, source of 'true' randomness for initializing random seed
    static std::random_device rd; 

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd()); 

    double sample;
    std::exponential_distribution<double> d(lam); 
    // get random number with normal distribution using gen as random source
    sample = d(gen);
	return sample;
}

// [[Rcpp::export]]
double myrunif(
	){
	// random device class instance, source of 'true' randomness for initializing random seed
    static std::random_device rd; 

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd()); 

    double sample;
    std::uniform_real_distribution<double> d(0.0, 1.0); 
    // get random number with normal distribution using gen as random source
    sample = d(gen);
	return sample;
}

// // [[Rcpp::export]]
// arma::cube to_cube(arma::vec values, int x, int y, int z) {
//   arma::vec v(values);
//   // v.randn();
//   
//   arma::cube res((const double*)v.begin(), x, y, z);
//   return res;
// }
// 
// // [[Rcpp::export]]
// arma::rowvec to_vec(arma::cube Q){
//   // arma::rowvec A = Q(arma::span(0), arma::span(1), arma::span::all);
//   arma::vec A = vectorise(Q);
//   return A.as_row();
// }


// [[Rcpp::export]]
arma::vec genpois(
	double tcurr,
	double T,
	double lam
	){
	double dt;
	int l = 0;
	int size = 2;
	arma::vec tg(size);
	while (tcurr <= T){
		dt = myrexp(lam);
		tcurr += dt;
		if(tcurr <= T){
			if (l >= size) {
				size = 2 * size;
				tg.resize(size); // amortization
			}
			tg(l) = tcurr;
			l++;
			//tg.print();
		}			 
	}
	tg = tg.head(l);
	return tg;
}


// [[Rcpp::export]]
arma::vec sampleChild(
	double tstart,
	double T,
	double lam,
	double b
	){
	double dt, bar, u;
	int size = 2;
	arma::vec child(size);
	int l = 0;
	double tcurr = tstart;
	double Tmax = 10.0 / lam + tstart;
	T = min(Tmax, T);
	while (tcurr <= T){
		dt = myrexp(lam*b);
		tcurr += dt;
		bar = exp(-lam * (tcurr - tstart));
		u = myrunif();
		if(tcurr <= T && u < bar){
			if (l >= size) {
				size = 2 * size;
				child.resize(size);
			}
			child(l) = tcurr;
			l++;
		}			 
	}
	child = child.head(l);
	return child;
}


// generate pair of nodes 

// [[Rcpp::export]]
arma::vec sampleHak(
	double T,
	double mu,
	double b,
	double lam
	){
	arma::vec parents = genpois(0, T, mu);
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

// [[Rcpp::export]]
arma::mat sampleBlockHak(
	double T,
	Rcpp::List A,
	arma::vec Z,
	arma::mat Mu,
	arma::mat B,
	double lam
	){
	int count = 0, size = 2;
	arma::mat alltimes(size, 3);
	int m = A.size();
	//printf("m: %d", m);
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
			temp = sampleHak(T, mu, b, lam);
			n_temp = temp.n_elem;
			// printf("n_temp: %d \n", n_temp);
			if(count + n_temp > size){
				size = max(2*size, count + n_temp);
				// size = count + n_temp;
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

// -----------------------------
// non-homogenuous baseline case

// [[Rcpp::export]]
arma::vec genpois_nonhomo(
	double t_start,
	double T,
	double window,
	arma::vec avec
	){
	int H = avec.n_elem, h;
	double period = window * H;
	double u, dt, tcurr = t_start, amax = avec.max();
	int l = 0, size = 2;
	arma::vec tg(size);
	while (tcurr <= T){
		dt = myrexp(amax);
		tcurr += dt;
		if(tcurr <= T){
			if (l >= size) {
				size = 2 * size;
				tg.resize(size); // amortization
			}
			h = floor((tcurr - floor(tcurr/period) * period)/window);// which window
			u = myrunif();
			if (u <= (avec(h)/amax)){
				tg(l) = tcurr;
				l++;
			}
			//tg.print();
		}			 
	}
	tg = tg.head(l);
	return tg;
}


// [[Rcpp::export]]
arma::vec sampleHak_nonhomo(
	double T,
	arma::vec avec,
	double b,
	double window,
	double lam
	){
	arma::vec parents = genpois_nonhomo(0, T, window, avec); // only change parents
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
arma::mat sampleBlockHak_nonhomo(
	double T,
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
	//printf("H: %d", H);
	//printf("m: %d", m);
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
			temp = sampleHak_nonhomo(T, avec, b, window, lam);
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
arma::mat sampleCCRMHak(
	double T,
	Rcpp::List A,
	arma::mat W1,
	arma::mat W2,
	double b,
	double lam
	){
	int count = 0, size = 2;
	int K = W1.n_cols;
	arma::mat alltimes(size, 3);
	int m = A.size();
	printf("m: %d", m);
	int n_edge, n_temp;
	double mu;
	int i,j,k,p,l;
	arma::vec temp;
	for (i = 0; i < m; i++){
		//printf("count: %d", count);
		arma::rowvec edge = A[i];
		n_edge = edge.n_elem;
		for (p = 0; p < n_edge; p++) {
			j = edge(p);
			mu = 0.0;
			for (l = 0; l < K; l++)
				mu += W1(i,l) * W2(j,l);

			temp = sampleHak(T, mu, b, lam);
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