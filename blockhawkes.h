#include "onlineblock.h"

double integral(double t, double T, double lam){
	return 1 - exp( - lam * (T - t));
}

double integral2(double t, double T1, double T2, double lam){
	return exp(- lam * (T1 - t)) - exp(- lam * (T2 - t));
}

double trigger(double t, double T, double lam){
	return lam * exp( - lam*(T - t));
}


Rcpp::List update_on(
	arma::mat tau,
	arma::mat Mu,
	arma::mat B,
	arma::rowvec Pi,
	arma::mat S,
	unordered_map<string, arma::vec> datamap,
	double t_start,
	double t_end,
	int m,
	int K,
	Rcpp::List A,
	double lam,
	double eta
	){
	arma::mat P1_mu(K,K), P2_mu(K,K), P1_B(K,K), P2_B(K,K), S_tp(m,K);
	arma::mat P1_mu_tp(K,K), P2_mu_tp(K,K), P1_B_tp(K,K), P2_B_tp(K,K), P_S_tp(m,K), Lambda(K,K);
	P1_mu.fill(0.0), P2_mu.fill(0.0), P1_B.fill(0.0), P2_B.fill(0.0), S_tp.fill(0.0), Lambda.fill(0.0);

	int l, k, n_edge;
	for (int i = 0; i < m; i++) {
		arma::rowvec edge = A[i];
		n_edge = edge.n_elem;
		for (int p = 0; p < n_edge; p++){
			int j = (int) edge(p);
			for (k = 0; k < K; k++) {
				for (l = 0; l < K; l++) {
					P2_mu(k,l) = P2_mu(k,l) + tau(i,k) * tau(j,l) * (t_end - t_start);
				}
			}
		}
	}

	unordered_map<string, arma::vec>:: iterator itr; 
	arma::vec timevec;
	string key;
	int ln,n;
    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
        key = itr->first;
        timevec = itr->second;
        arma::vec index = split(key);
        int i = (int) index(0), j = (int) index(1);

        P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0), P1_B_tp.fill(0.0), P2_B_tp.fill(0.0), P_S_tp.fill(0.0);
        P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
        	double t_current = timevec(n);
        	if (t_current > t_start) {
        		double intensity = 0;
        		Lambda.fill(0.0); // store temporary intensity values
        		for (int n1 = 0; n1 < n; n1++) {
        			double t1 = timevec(n1);
        			intensity += trigger(t1, t_current, lam);
        		}
        		P2_B_tp = P2_B_tp + integral(t_current, t_end, lam);
        		for (k = 0; k < K; k++){
        			for (l = 0; l < K; l++){
        				Lambda(k,l) = Mu(k,l) + B(k,l) * intensity;
        			}
        		}
        		for (k = 0; k < K; k++) {
        			for (l = 0; l < K; l++){
        				P1_mu_tp(k,l) += 1.0 / Lambda(k,l);
        				P1_B_tp(k,l) += intensity / Lambda(k,l);
        				P_S_tp(k,l) += log(Lambda(k,l));
        			}
        		}
        	} else {
        		P2_B_tp = P2_B_tp + integral2(t_current, t_start, t_end, lam);
        	}
        }

        for (k = 0; k < K; k++) {
        	for (l = 0; l < K; l++) {
        		P1_mu(k,l) += tau(i,k) * tau(j,l) * P1_mu_tp(k,l);
        		P1_B(k,l) += tau(i,k) * tau(j,l) * P1_B_tp(k,l);
        		P2_B(k,l) += tau(i,k) * tau(j,l) * P2_B_tp(k,l);
        	}
        }

        // update S
        for (k = 0; k < K; k++) {
        	for (l = 0; l < K; l++) {
        		S_tp(i,k) += tau(j,l) * (P_S_tp(k,l) - B(k,l) * P2_B_tp(k,l));
        	}
        }
    } 

    // update S, second part
    for (int i = 0; i < m; i++) {
    	arma::rowvec edge = A[i];
    	n_edge = edge.n_elem;
    	for (k = 0; k < K; k++) {
    		for (int p = 0; p < n_edge; p++) {
    			int j = (int) edge(p);
    			for (l = 0; l < K; l++) {
    				S_tp(i,k) = S_tp(i,k) - tau(j,l) * Mu(k,l) * (t_end - t_start);
    			}
    		}
    	}
    }


    // update parameters
    S = S + S_tp;
    arma::mat grad_B = P1_B - P2_B;
    //grad_B.print();
    arma::mat grad_mu = P1_mu - P2_mu;
    //grad_mu.print();
    arma::mat B_new = B + eta * grad_B;
    //printf("B new is: \n");   
    //B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //printf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
    	for (l = 0; l < K; l++) {
    		if (B_new(k,l) <= 0.0) 
    			B_new(k,l) = B(k,l) / 2.0;
    		if (B_new(k,l) > 2 * B(k,l))
    			B_new(k,l) = B(k,l) * 2.0;
    		if (Mu_new(k,l) <= 0.0)
    			Mu_new(k,l) = Mu(k,l) / 2.0;
    		if (Mu_new(k,l) > 2 * Mu(k,l))
    			Mu_new(k,l) = Mu(k,l) * 2.0;
    	}
    }

    arma::mat tau_new(m,K);
    tau_new.fill(0.0);

    arma::rowvec s;
    for (int i = 0; i < m; i++) {
    	s = S.row(i) - S.row(i).max();
    	s = s + log(Pi + 0.000001);
    	s = exp(s)/sum(exp(s));
    	tau_new.row(i) = s;
    }

    for (k = 0; k < K; k++) {
    	Pi(k) = sum(tau.col(k)) / (m + 0.0);
    }	

	return Rcpp::List::create(Rcpp::Named("tau") = tau_new,
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("S") = S);
}


// [[Rcpp::export]]
Rcpp::List online_estimator(
	arma::mat alltimes,
	Rcpp::List A,
	int m,
	int K,
	double T,
	double dT,
	double lam,
	arma::mat B_start,
	arma::mat Mu_start,
	arma::mat tau_start
	){
	// initialization
	arma::rowvec Pi(K);
	Pi.fill(1.0 / K);
	arma::mat B(K,K), Mu(K,K), S(m,K);
	arma::mat tau(m,K);
	B.fill(0.5), Mu.fill(0.5), S.fill(0.0);
	//B = B_start, Mu = Mu_start;
	for (int i = 0; i < m; i++) {
		arma::rowvec tt(K);
		for (int k = 0; k < K; k++) {
			tt(k) = myrunif();
		}
		tt = tt / sum(tt);
		tau.row(i) = tt;
	}
	//tau = tau_start;

	int nall = alltimes.n_rows;
	int trunc_pos = 0, start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
	int N = floor(T / dT);
	int nsave = floor(5.0 / dT);
	queue<double> trunc_pos_queue;

	double Tn, t_current, t_start, eta;
	arma::rowvec event; 
	arma::mat truncdata;
	Rcpp::List paralist;
	unordered_map<string, arma::vec> datamap;
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
		if ( (int) trunc_pos_queue.size() == nsave) {
			trunc_pos = trunc_pos_queue.front();
			trunc_pos_queue.pop();
		}
		truncdata = alltimes.rows(trunc_pos, end_pos - 1);
		datamap = transfer(truncdata);
		t_start = Tn - dT;
		ln_curr = end_pos;
		n_t = ln_curr - ln_prev;
		eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
		paralist = update_on(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
		arma::mat tau_new = paralist["tau"], Mu_new = paralist["Mu"], B_new = paralist["B"], S_new = paralist["S"];
		arma::rowvec Pi_new = paralist["Pi"];
		tau = tau_new; 
		Mu = Mu_new, B = B_new, S = S_new, Pi = Pi_new;
		trunc_pos_queue.push(start_pos);
  		start_pos = curr_pos;
  		ln_prev = ln_curr;
  		printf("iter: %d; number: %d \n", n, n_t); 
  		B.print();
  		Mu.print();
  		//S.print();
  		printf("=============\n");
	}

	return Rcpp::List::create(
                          Rcpp::Named("Mu") = Mu,
                          Rcpp::Named("B") = B,
                          Rcpp::Named("Pi") = Pi);
}