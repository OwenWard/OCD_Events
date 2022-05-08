#include "onlineblock.h"

double integral(double t, double T, double lam){
	return 1 - exp( - lam * (T - t));
}

double integral2(double t, double T1, double T2, double lam){
	return exp(- lam * (T1 - t)) - exp(- lam * (T2 - t));
}

double trigger(double t, double T, double lam){
	return lam * exp( - lam * (T - t));
}

double trigger_lam(double t, double T, double lam){
    return exp( - lam * (T - t));
}

// arma::rowvec correct_tau(arma::rowvec tau){
//   arma::rowvec mod_tau = tau;
//   int m = tau.n_cols;
//   for(int j = 0; j < m; ++j) {
//     mod_tau[j] = min(tau[j], 1-1e-7);
//     mod_tau[j] = max(tau[j], 1e-7);
//   }
//   // mod_tau = max(tau, 1e-7);
//   // then normalise to sum to 1 again
//   mod_tau = mod_tau/sum(mod_tau);
//   return mod_tau;
// }


double ELBO_Hak(
    unordered_map<string, arma::vec> datamap,
    double t_start,
    double t_end,
    arma::mat tau,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    Rcpp::List A,
    double lam,
    int m,
    int K
    ){
    double elbo = 0;

    int l, k, n_edge;
    arma::vec timevec;
    string key;
    int ln,n;
    double intensity;
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
            if (t_current > t_start) {
                intensity = 0.0;
                Lambda.fill(0.0); // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                }
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) = Mu(k,l) + B(k,l) * intensity;
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

    // mu * T
    for (int i = 0; i < m; i++) {
        arma::rowvec edge = A[i];
        n_edge = edge.n_elem;
        for (int p = 0; p < n_edge; p++) {
            int j = (int) edge(p);
            for (k = 0; k < K; k++) {
                for (l = 0; l < K; l++) {
                    elbo -= tau(i,k) * tau(j,l) * Mu(k,l)* (t_end - t_start);
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
double get_elbo_Hak(
    arma::mat alltimes,
    double t_start,
    double t_end,
    arma::mat tau,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    Rcpp::List A,
    double lam,
    int m,
    int K
    ){
    unordered_map<string, arma::vec> datamap = transfer(alltimes);
    double elbo = ELBO_Hak(datamap, t_start, t_end, tau, Mu, B, Pi, A, lam, m, K);
    return elbo;
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

        if (i == j)
            continue;

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
    //Rprintf("B new is: \n");   
    //B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //Rprintf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
    	for (l = 0; l < K; l++) {
    		if (B_new(k,l) <= 0.0) {
    			B_new(k,l) = B(k,l) / 2.0;
            } else if (B_new(k,l) > 2 * B(k,l)) {
                B_new(k,l) = B(k,l) * 2.0;
            }
    		
    		if (Mu_new(k,l) <= 0.0) {
    			Mu_new(k,l) = Mu(k,l) / 2.0;
            } else if (Mu_new(k,l) > 2 * Mu(k,l)) {
                Mu_new(k,l) = Mu(k,l) * 2.0;
            }

    	}
    }

    arma::mat tau_new(m,K);
    tau_new.fill(0.0);
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
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("S") = S);
}


Rcpp::List update_lam_trunc(
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
    double eta,
    double trunc_length
    ){
    arma::mat P1_mu(K,K), P2_mu(K,K), P1_B(K,K), P2_B(K,K), S_tp(m,K);
    arma::mat P1_mu_tp(K,K), P2_mu_tp(K,K), P1_B_tp(K,K), P2_B_tp(K,K);
    arma::mat P_S_tp(m,K), Lambda(K,K);
    P1_mu.fill(0.0), P2_mu.fill(0.0), P1_B.fill(0.0);
    P2_B.fill(0.0), S_tp.fill(0.0), Lambda.fill(0.0);
    arma::mat lam_store(K,K);

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
    double R = trunc_length / lam;

    double intensity_lam1, intensity_lam2, intensity, grad_lam = 0.0;
    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stores the value part 
        
        key = itr->first;
        timevec = itr->second;
        arma::vec index = split(key);
        int i = (int) index(0), j = (int) index(1);

        if (i == j)
            continue;
        
        P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0);
        P1_B_tp.fill(0.0), P2_B_tp.fill(0.0), P_S_tp.fill(0.0);
        lam_store.fill(0.0);
        P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = 0, intensity_lam1 = 0.0, intensity_lam2 = 0.0;
                Lambda.fill(0.0); // store temporary intensity values
                for (int n1 = n - 1; n1 >= 0; n1--) {
                    double t1 = timevec(n1);
                    if (t_current - t1 > R)
                        break;
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_B_tp = P2_B_tp + integral(t_current, t_end, lam);
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) = Mu(k,l) + B(k,l) * intensity;
                    }
                }
                lam_store = lam_store + B * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - B * ((t_end - t_current) * exp(-lam*(t_end - t_current)));
                
                for (k = 0; k < K; k++) {
                    for (l = 0; l < K; l++){
                        P1_mu_tp(k,l) += 1.0 / Lambda(k,l);
                        P1_B_tp(k,l) += intensity / Lambda(k,l);
                        P_S_tp(k,l) += log(Lambda(k,l));
                    }
                }
            } else {
                P2_B_tp = P2_B_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + B * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
                
            }
        }
        // Rprintf("lam store \n");
        // sum(lam_store, 0).print();
        // lam_store.fill(1.0);
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
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
        // cout<<grad_lam<<endl;
        // Rprintf("Check S_tp here \n");
        // sum(S_tp, 0).print();
        // // Rprintf("Check tau \n");
        // // sum(tau, 0).print();
        // // Rprintf("Check P_S_tp \n");
        // // sum(P_S_tp, 0).print();
        // Rprintf("Check P2_B_tp \n");
        // sum(P2_B_tp, 0).print();
    } 
    
    
    // Rprintf("Check B \n");
    // B.print();
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
    // Rprintf("P1_B \n");
    // P1_B.print();
    // Rprintf("P2_B \n");
    // P2_B.print();
    // Rprintf("------- \n");
    // grad_B.print();
    arma::mat grad_mu = P1_mu - P2_mu;
    //grad_mu.print();
    arma::mat B_new = B + eta * grad_B;
    // Rprintf("B new is: \n");
    // B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //Rprintf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
            if (B_new(k,l) <= 0.0) 
                B_new(k,l) = B(k,l) / 2.0;
            else if (B_new(k,l) > 2 * B(k,l))
                B_new(k,l) = B(k,l) * 2.0;
            if (Mu_new(k,l) <= 0.0)
                Mu_new(k,l) = Mu(k,l) / 2.0;
            else if (Mu_new(k,l) > 2 * Mu(k,l))
                Mu_new(k,l) = Mu(k,l) * 2.0;
        }
    }
    // Rprintf("Grad lambda \n");
    // cout << grad_lam << endl;
    // double lam_new = lam + eta * grad_lam;
    // if (lam_new > 1.5*lam) {
    //     lam_new = 1.5 * lam;
    // } else if (lam_new <= 0.0) {
    //     lam_new = lam/2.0;
    // }
    double lam_new = 0.15;

    arma::mat tau_new(m,K);
    
    tau_new.fill(0.0);
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
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam_new,
                          Rcpp::Named("S") = S);
}

// declaration
Rcpp::List update_lam(
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
    );


// stochastic version

Rcpp::List update_lam_stoch(
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
    double eta,
    double percent
    ){
    arma::mat P1_mu(K,K), P2_mu(K,K), P1_B(K,K), P2_B(K,K), S_tp(m,K);
    arma::mat P1_mu_tp(K,K), P2_mu_tp(K,K), P1_B_tp(K,K), P2_B_tp(K,K), P_S_tp(m,K), Lambda(K,K);
    P1_mu.fill(0.0), P2_mu.fill(0.0), P1_B.fill(0.0), P2_B.fill(0.0), S_tp.fill(0.0), Lambda.fill(0.0);
    arma::mat lam_store(K,K);


    int l, k;

    unordered_map<string, arma::vec>:: iterator itr; 
    arma::vec timevec;
    string key;
    int ln,n;
    double intensity_lam1, intensity_lam2, intensity, grad_lam = 0.0;
    double u_var;

    int zero_count = 0;

    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
        u_var = myrunif();
        if (u_var > percent) {
            continue;
        }

        key = itr->first;
        timevec = itr->second;
        arma::vec index = split(key);
        int i = (int) index(0), j = (int) index(1);

        if (i == j)
            continue;
        
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                P2_mu(k,l) = P2_mu(k,l) + tau(i,k) * tau(j,l) * (t_end - t_start);
            }
        }

        // update S, second part
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                S_tp(i,k) = S_tp(i,k) - tau(j,l) * Mu(k,l) * (t_end - t_start);
            }
            
        }


        P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0), P1_B_tp.fill(0.0), P2_B_tp.fill(0.0), P_S_tp.fill(0.0);
        lam_store.fill(0.0);
        P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        if (ln == 0) {
            zero_count ++;
        }
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda.fill(eps); // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_B_tp = P2_B_tp + integral(t_current, t_end, lam);
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) += Mu(k,l) + B(k,l) * intensity;
                    }
                }
                lam_store = lam_store + B * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - B * ((t_end - t_current) * exp(-lam*(t_end - t_current)));
                for (k = 0; k < K; k++) {
                    for (l = 0; l < K; l++){
                        P1_mu_tp(k,l) += 1.0 / Lambda(k,l);
                        P1_B_tp(k,l) += intensity / Lambda(k,l);
                        P_S_tp(k,l) += log(Lambda(k,l));
                    }
                }
            } else {
                P2_B_tp = P2_B_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + B * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
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

    Rprintf("zero_count: %d \n", zero_count);
    cout<<"Print this?"<<endl;


    // update parameters
    S = S + S_tp;
    arma::mat grad_B = P1_B - P2_B;
    //grad_B.print();
    arma::mat grad_mu = P1_mu - P2_mu;
    //grad_mu.print();
    arma::mat B_new = B + eta * grad_B;
    //Rprintf("B new is: \n");   
    //B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //Rprintf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
            if (B_new(k,l) <= 0.0) 
                B_new(k,l) = B(k,l) / 2.0;
            else if (B_new(k,l) > 2 * B(k,l))
                B_new(k,l) = B(k,l) * 2.0;
            if (Mu_new(k,l) <= 0.0)
                Mu_new(k,l) = Mu(k,l) / 2.0;
            else if (Mu_new(k,l) > 2 * Mu(k,l))
                Mu_new(k,l) = Mu(k,l) * 2.0;
        }
    }
    double lam_new = lam + eta * grad_lam;
    if (lam_new > 5*lam) {
        lam_new = 5 * lam;
    } else if (lam_new <= 0.0) {
        lam_new = lam/2.0;
    }

    arma::mat tau_new(m,K);
    tau_new.fill(0.0);

    for (int i = 0; i < m; i++) {
      arma::rowvec s = arma::log(Pi) + S.row(i);
      s = s - max(s);
      s = exp(s)/sum(exp(s));
      tau_new.row(i) = correct_tau(s);
    }

    for (k = 0; k < K; k++) {
        Pi(k) = sum(tau_new.col(k)) / (m + eps);
    }   

    return Rcpp::List::create(Rcpp::Named("tau") = tau_new,
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam_new,
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
	arma::mat tau_start,
    double percent = 1.0,
    bool is_elbo = false
	){
	// initialization
	arma::rowvec Pi(K);
	Pi.fill(1.0 / K);
	arma::mat B(K,K), Mu(K,K), S(m,K);
	arma::mat tau(m,K);
	tau.fill(1.0/K);
	S.fill(1.0/K);
	for (int k = 0; k < K; k++) {
        for (int l=0; l < K; l++) {
            B(k,l) = myrunif();
            Mu(k,l) = myrunif();
        }
    }
    //B.fill(0.5), Mu.fill(0.5); 
    // S.fill(0.0);
	//B = B_start, Mu = Mu_start;
	// for (int i = 0; i < m; i++) {
	// 	arma::rowvec tt(K);
	// 	for (int k = 0; k < K; k++) {
	// 		tt(k) = myrunif();
	// 	}
	// 	tt = tt / sum(tt);
	// 	tau.row(i) = tt;
	// }
	//tau = tau_start;

	int nall = alltimes.n_rows;
	int trunc_pos = 0, start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
	int N = floor(T / dT);
	int nsave = floor(5.0 / dT);
	queue<double> trunc_pos_queue;
    
    arma::vec elbo_vec(N);
    double elbo = 0;
    arma::mat prevdata;

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
		// datamap = transfer(truncdata);
        datamap = transfer2(truncdata, A, m);
		t_start = Tn - dT;
		ln_curr = end_pos;
		n_t = ln_curr - ln_prev;
		eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K) / percent;
        paralist = update_lam(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
        // paralist = update_lam_stoch(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta, percent);
		arma::mat tau_new = paralist["tau"], Mu_new = paralist["Mu"], B_new = paralist["B"], S_new = paralist["S"];
		arma::rowvec Pi_new = paralist["Pi"];
        double lam_new = paralist["lam"];
		tau = tau_new; 
		Mu = Mu_new, B = B_new, S = S_new, Pi = Pi_new;
        lam = lam_new;
		trunc_pos_queue.push(start_pos);
  		start_pos = curr_pos;
  		ln_prev = ln_curr;
  		Rprintf("iter: %d; number: %d \n", n, n_t); 
  		B.print();
  		Mu.print();
        Rprintf("lam: %2.3f", lam);

        if (is_elbo) {
            prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
            elbo = get_elbo_Hak(prevdata, 0.0, Tn, tau, Mu, B, Pi, A, lam, m, K);
            elbo_vec(n) = elbo / ln_curr;
        }

  		//S.print();
  		Rprintf("=============\n");
	}

	return Rcpp::List::create(
                          Rcpp::Named("Mu") = Mu,
                          Rcpp::Named("B") = B,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam,
                          Rcpp::Named("tau") = tau,
                          Rcpp::Named("elbo") = elbo_vec);
}



Rcpp::List update_lam(
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
    arma::mat lam_store(K,K);


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
    double intensity_lam1, intensity_lam2, intensity, grad_lam = 0.0;
    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
        key = itr->first;
        timevec = itr->second;
        arma::vec index = split(key);
        int i = (int) index(0), j = (int) index(1);

        if (i == j)
            continue;
        
        P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0), P1_B_tp.fill(0.0), P2_B_tp.fill(0.0), P_S_tp.fill(0.0);
        lam_store.fill(0.0);
        P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda.fill(eps); // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_B_tp = P2_B_tp + integral(t_current, t_end, lam);
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) += Mu(k,l) + B(k,l) * intensity;
                    }
                }
                lam_store = lam_store + B * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - B * ((t_end - t_current) * exp(-lam*(t_end - t_current)));
                for (k = 0; k < K; k++) {
                    for (l = 0; l < K; l++){
                        P1_mu_tp(k,l) += 1.0 / Lambda(k,l);
                        P1_B_tp(k,l) += intensity / Lambda(k,l);
                        P_S_tp(k,l) += log(Lambda(k,l));
                    }
                }
            } else {
                P2_B_tp = P2_B_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + B * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
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
    //Rprintf("B new is: \n");   
    //B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //Rprintf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
            if (B_new(k,l) <= 0.0) 
                B_new(k,l) = B(k,l) / 2.0;
            else if (B_new(k,l) > 2 * B(k,l))
                B_new(k,l) = B(k,l) * 2.0;
            if (Mu_new(k,l) <= 0.0)
                Mu_new(k,l) = Mu(k,l) / 2.0;
            else if (Mu_new(k,l) > 2 * Mu(k,l))
                Mu_new(k,l) = Mu(k,l) * 2.0;
        }
    }
    double lam_new = lam + eta * grad_lam;
    if (lam_new > 5*lam) {
        lam_new = 5 * lam;
    } else if (lam_new <= 0.0) {
        lam_new = lam/2.0;
    }

    arma::mat tau_new(m,K);
    tau_new.fill(0.0);

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
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam_new,
                          Rcpp::Named("S") = S);
}


double get_grad_lam(
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
    double lam
    ){
    arma::mat lam_store(K,K), Lambda(K,K);
    unordered_map<string, arma::vec>:: iterator itr; 
    arma::vec timevec;
    string key;
    int ln, n, k, l;
    double intensity_lam1, intensity_lam2, intensity, grad_lam = 0.0;

    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
        key = itr->first;
        timevec = itr->second;
        arma::vec index = split(key);
        int i = (int) index(0), j = (int) index(1);

        lam_store.fill(0.0);

        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = 0, intensity_lam1 = 0.0, intensity_lam2 = 0.0;
                Lambda.fill(0.0); // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) = Mu(k,l) + B(k,l) * intensity;
                    }
                }
                lam_store = lam_store + B * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - B * ((t_end - t_current) * exp(-lam*(t_end - t_current)));
            } else {
                lam_store = lam_store + B * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
            }
        }        

    }
    return grad_lam;
}

// [[Rcpp::export]]
double test_lam(
    arma::mat tau,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    arma::mat S,
    arma::mat alltimes,
    double t_start,
    double t_end,
    int m,
    int K,
    Rcpp::List A,
    double lam
    ){
    unordered_map<string, arma::vec> datamap = transfer(alltimes);
    double grad = get_grad_lam(tau, Mu, B, Pi, S, datamap, t_start, t_end, m, K, A, lam);
    return grad;
}




Rcpp::List update_lam_eff(
    arma::mat tau,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    arma::mat S,
    unordered_map<string, std::deque<double>> &datamap,
    double t_start,
    double t_end,
    int m,
    int K,
    //Rcpp::List A,
    double lam,
    double eta
    ){
    arma::mat P1_mu(K,K), P2_mu(K,K), P1_B(K,K), P2_B(K,K), S_tp(m,K);
    arma::mat P1_mu_tp(K,K), P2_mu_tp(K,K), P1_B_tp(K,K), P2_B_tp(K,K), P_S_tp(m,K), Lambda(K,K);
    P1_mu.fill(0.0), P2_mu.fill(0.0), P1_B.fill(0.0), P2_B.fill(0.0), S_tp.fill(0.0), Lambda.fill(0.0);
    arma::mat lam_store(K,K);


    int l, k;
    //int n_edge;

    /*
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
    */

    unordered_map<string, std::deque<double>>:: iterator itr; 
    arma::vec timevec;
    std::deque<double> timeque;
    string key;
    int ln,n;
    double intensity_lam1, intensity_lam2, intensity, grad_lam = 0.0;
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
                P2_mu(k,l) = P2_mu(k,l) + tau(i,k) * tau(j,l) * (t_end - t_start);
            }
        }

        // update S, second part
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                S_tp(i,k) = S_tp(i,k) - tau(j,l) * Mu(k,l) * (t_end - t_start);
            }
            
        }

        P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0), P1_B_tp.fill(0.0), P2_B_tp.fill(0.0), P_S_tp.fill(0.0);
        lam_store.fill(0.0);
        //P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda.fill(eps); // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_B_tp = P2_B_tp + integral(t_current, t_end, lam);
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) += Mu(k,l) + B(k,l) * intensity;
                    }
                }
                lam_store = lam_store + B * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - B * ((t_end - t_current) * exp(-lam*(t_end - t_current)));
                for (k = 0; k < K; k++) {
                    for (l = 0; l < K; l++){
                        P1_mu_tp(k,l) += 1.0 / Lambda(k,l);
                        P1_B_tp(k,l) += intensity / Lambda(k,l);
                        P_S_tp(k,l) += log(Lambda(k,l));
                    }
                }
            } else {
                P2_B_tp = P2_B_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + B * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
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
    /*
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
    */

    // update parameters
    S = S + S_tp;
    arma::mat grad_B = P1_B - P2_B;
    //grad_B.print();
    arma::mat grad_mu = P1_mu - P2_mu;
    //grad_mu.print();
    arma::mat B_new = B + eta * grad_B;
    //Rprintf("B new is: \n");   
    //B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //Rprintf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
            if (B_new(k,l) <= 0.0) 
                B_new(k,l) = B(k,l) / 2.0;
            else if (B_new(k,l) > 2 * B(k,l))
                B_new(k,l) = B(k,l) * 2.0;
            if (Mu_new(k,l) <= 0.0)
                Mu_new(k,l) = Mu(k,l) / 2.0;
            else if (Mu_new(k,l) > 2 * Mu(k,l))
                Mu_new(k,l) = Mu(k,l) * 2.0;
        }
    }
    double lam_new = lam + eta * grad_lam;
    if (lam_new > 5*lam) {
        lam_new = 5 * lam;
    } else if (lam_new <= 0.0) {
        lam_new = lam/2.0;
    }

    arma::mat tau_new(m,K);
    tau_new.fill(0.0);

    for (int i = 0; i < m; i++) {
      arma::rowvec s = arma::log(Pi) + S.row(i);
      s = s - max(s);
      s = exp(s)/sum(exp(s));
      tau_new.row(i) = correct_tau(s);
    }

    for (k = 0; k < K; k++) {
        Pi(k) = sum(tau.col(k)) / (m + eps);
    }   

    return Rcpp::List::create(Rcpp::Named("tau") = tau_new,
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam_new,
                          Rcpp::Named("S") = S);
}



Rcpp::List update_lam_eff_revised(
    arma::mat tau,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    arma::mat S,
    unordered_map<string, std::deque<double>> &datamap,
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
    arma::mat lam_store(K,K);


    int l, k;
    
    int n_edge;

    for (int i = 0; i < m; i++) {
        arma::rowvec edge = A[i];
        n_edge = edge.n_elem;
        for (int p = 0; p < n_edge; p++){
            int j = (int) edge(p);
            if(i != j){
              for (k = 0; k < K; k++) {
                for (l = 0; l < K; l++) {
                  P2_mu(k,l) = P2_mu(k,l) + tau(i,k) * tau(j,l) * (t_end - t_start);
                }
              }
              
              // update S, second part
              for (k = 0; k < K; k++) {
                for (l = 0; l < K; l++) {
                  S_tp(i,k) = S_tp(i,k) - tau(j,l) * Mu(k,l) * (t_end - t_start);
                }
                
              }
            }
            
        }
    }


    unordered_map<string, std::deque<double>>:: iterator itr; 
    arma::vec timevec;
    std::deque<double> timeque;
    string key;
    int ln,n;
    double intensity_lam1, intensity_lam2, intensity, grad_lam = 0.0;
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
                P2_mu(k,l) = P2_mu(k,l) + tau(i,k) * tau(j,l) * (t_end - t_start);
            }
        }

        // update S, second part
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                S_tp(i,k) = S_tp(i,k) - tau(j,l) * Mu(k,l) * (t_end - t_start);
            }
            
        }
        */

        P1_mu_tp.fill(0.0), P2_mu_tp.fill(0.0), P1_B_tp.fill(0.0), P2_B_tp.fill(0.0), P_S_tp.fill(0.0);
        lam_store.fill(0.0);
        //P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda.fill(eps); // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_B_tp = P2_B_tp + integral(t_current, t_end, lam);
                for (k = 0; k < K; k++){
                    for (l = 0; l < K; l++){
                        Lambda(k,l) += Mu(k,l) + B(k,l) * intensity;
                    }
                }
                lam_store = lam_store + B * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - B * ((t_end - t_current) * exp(-lam*(t_end - t_current)));
                for (k = 0; k < K; k++) {
                    for (l = 0; l < K; l++){
                        P1_mu_tp(k,l) += 1.0 / Lambda(k,l);
                        P1_B_tp(k,l) += intensity / Lambda(k,l);
                        P_S_tp(k,l) += log(Lambda(k,l));
                    }
                }
            } else {
                P2_B_tp = P2_B_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + B * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
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

    /*
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
    */

    // update parameters
    S = S + S_tp;
    arma::mat grad_B = P1_B - P2_B;
    //grad_B.print();
    arma::mat grad_mu = P1_mu - P2_mu;
    //grad_mu.print();
    arma::mat B_new = B + eta * grad_B;
    //Rprintf("B new is: \n");   
    //B_new.print();
    arma::mat Mu_new = Mu + eta * grad_mu;
    //Rprintf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    for (k = 0; k < K; k++) {
        for (l = 0; l < K; l++) {
            if (B_new(k,l) <= 0.0) 
                B_new(k,l) = B(k,l) / 2.0;
            else if (B_new(k,l) > 2 * B(k,l))
                B_new(k,l) = B(k,l) * 2.0;
            if (Mu_new(k,l) <= 0.0)
                Mu_new(k,l) = Mu(k,l) / 2.0;
            else if (Mu_new(k,l) > 2 * Mu(k,l))
                Mu_new(k,l) = Mu(k,l) * 2.0;
        }
    }
    double lam_new = lam + eta * grad_lam;
    if (lam_new > 5*lam) {
        lam_new = 5 * lam;
    } else if (lam_new <= 0.0) {
        lam_new = lam/2.0;
    }

    arma::mat tau_new(m,K);
    tau_new.fill(0.0);

    for (int i = 0; i < m; i++) {
      arma::rowvec s = arma::log(Pi) + S.row(i);
      s = s - max(s);
      s = exp(s)/sum(exp(s));
      tau_new.row(i) = correct_tau(s);
    }

    for (k = 0; k < K; k++) {
        Pi(k) = sum(tau_new.col(k)) / (m + eps);
    }   

    return Rcpp::List::create(Rcpp::Named("tau") = tau_new,
                          Rcpp::Named("Mu") = Mu_new,
                          Rcpp::Named("B") = B_new,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam_new,
                          Rcpp::Named("S") = S);
}


// [[Rcpp::export]]
Rcpp::List online_estimator_eff(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    double lam,
    arma::mat B_start,
    arma::mat Mu_start,
    arma::mat tau_start,
    bool is_elbo = false
    ){

    // create empty map
    unordered_map<string, std::deque<double>> datamap;
    datamap = transfer_create(A, m);

    // initialization
    arma::rowvec Pi(K);
    Pi.fill(1.0 / K);
    arma::mat B(K,K), Mu(K,K), S(m,K);
    arma::mat tau(m,K);
    for (int k = 0; k < K; k++) {
        for (int l=0; l < K; l++) {
            B(k,l) = myrunif();
            Mu(k,l) = myrunif();
        }
    }
    //B.fill(0.5), Mu.fill(0.5); 
    S.fill(0.0);
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
    int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
    int N = floor(T / dT);

    double R = 5.0;
    
    arma::vec elbo_vec(N);
    double elbo = 0;
    arma::mat prevdata;

    double Tn, t_current, t_start, eta;
    arma::rowvec event; 
    arma::mat truncdata;
    Rcpp::List paralist;

    for (int n = 0; n < N; n++ ){
        // R = min(5.0 / lam, 10.0);
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
       
        
       
        //datamap = transfer_eff2(datamap, truncdata, R);
        
        transfer_eff(datamap, truncdata, R);
        
        t_start = Tn - dT;
        ln_curr = end_pos;
        n_t = ln_curr - ln_prev;
        eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
        // paralist = update_lam_eff(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
        paralist = update_lam_eff(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, lam, eta);
        arma::mat tau_new = paralist["tau"], Mu_new = paralist["Mu"], B_new = paralist["B"], S_new = paralist["S"];
        arma::rowvec Pi_new = paralist["Pi"];
        double lam_new = paralist["lam"];
        tau = tau_new; 
        Mu = Mu_new, B = B_new, S = S_new, Pi = Pi_new;
        lam = lam_new;
        start_pos = curr_pos;
        ln_prev = ln_curr;
        Rprintf("iter: %d; number: %d \n", n, n_t); 
        B.print();
        Mu.print();
        Rprintf("lam: %2.3f", lam);

        if (is_elbo) {
            prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
            elbo = get_elbo_Hak(prevdata, 0.0, Tn, tau, Mu, B, Pi, A, lam, m, K);
            elbo_vec(n) = elbo / ln_curr;
        }

        //S.print();
        Rprintf("=============\n");
    }

    return Rcpp::List::create(
                          Rcpp::Named("Mu") = Mu,
                          Rcpp::Named("B") = B,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam,
                          Rcpp::Named("tau") = tau,
                          Rcpp::Named("elbo") = elbo_vec);
}


// [[Rcpp::export]]
Rcpp::List online_estimator_eff_revised(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    double lam,
    arma::mat B_start,
    arma::mat Mu_start,
    arma::mat tau_start,
    int inter_T,
    bool is_elbo = false
    ){

    // create empty map
    unordered_map<string, std::deque<double>> datamap;
    datamap = transfer_create_empty();

    // initialization
    arma::rowvec Pi(K);
    Pi.fill(1.0 / K);
    arma::mat B(K,K), Mu(K,K), S(m,K);
    arma::mat tau(m,K);
    for (int k = 0; k < K; k++) {
        for (int l=0; l < K; l++) {
            B(k,l) = myrunif();
            Mu(k,l) = myrunif();
        }
    }
    //B.fill(0.5), Mu.fill(0.5); 
    S.fill(1.0/K);
    tau.fill(1.0/K);
    //B = B_start, Mu = Mu_start;
    // for (int i = 0; i < m; i++) {
    //     arma::rowvec tt(K);
    //     for (int k = 0; k < K; k++) {
    //         tt(k) = myrunif();
    //     }
    //     tt = tt / sum(tt);
    //     tau.row(i) = tt;
    // }
    //tau = tau_start;

    int nall = alltimes.n_rows;
    int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
    int N = floor(T / dT);
    int slices = int(N/inter_T);
    int ind = 0;

    double R = 5.0;
    
    arma::vec elbo_vec(N);
    double elbo = 0;
    arma::mat prevdata;

    double Tn, t_current, t_start, eta;
    arma::rowvec event; 
    arma::mat truncdata;
    Rcpp::List paralist;
    
    // to inspect intermediate values of clustering
    arma::cube inter_tau(m,K,slices+1);
    arma::cube inter_B(K,K,N);
    arma::cube inter_mu(K,K,N);

    for (int n = 0; n < N; n++ ){
        // R = min(5.0 / lam, 10.0);
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

        //datamap = transfer_eff2(datamap, truncdata, R);
        // transfer_eff(datamap, truncdata, R);
        // Rprintf("Datamap \n");
        transfer_dynamic(datamap, truncdata, R, Tn);

        t_start = Tn - dT;
        ln_curr = end_pos;
        n_t = ln_curr - ln_prev;
        eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
        // paralist = update_lam_eff(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
        paralist = update_lam_eff_revised(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);

        arma::mat tau_new = paralist["tau"], Mu_new = paralist["Mu"], B_new = paralist["B"], S_new = paralist["S"];
        arma::rowvec Pi_new = paralist["Pi"];
        double lam_new = paralist["lam"];
        tau = tau_new; 
        Mu = Mu_new, B = B_new, S = S_new, Pi = Pi_new;
        lam = lam_new;
        start_pos = curr_pos;
        ln_prev = ln_curr;
        //Rprintf("iter: %d; number: %d \n", n, n_t); 
        inter_mu.slice(n) = Mu;
        inter_B.slice(n) = B;
        //B.print();
        //Mu.print();
        //Rprintf("lam: %2.3f", lam);
        
        
        if(n % inter_T == 0 ){
          inter_tau.slice(ind) = tau;
          ind = ind + 1;
        }

        if (is_elbo) {
            prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
            elbo = get_elbo_Hak(prevdata, 0.0, Tn, tau, Mu, B, Pi, A, lam, m, K);
            elbo_vec(n) = elbo / ln_curr;
        }

        //S.print();
        //Rprintf("=============\n");
    }

    return Rcpp::List::create(
                          Rcpp::Named("Mu") = Mu,
                          Rcpp::Named("B") = B,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam,
                          Rcpp::Named("tau") = tau,
                          Rcpp::Named("early_tau")= inter_tau,
                          Rcpp::Named("inter_B") = inter_B,
                          Rcpp::Named("inter_mu") = inter_mu,
                          Rcpp::Named("elbo") = elbo_vec);
}



// batch optimization
// [[Rcpp::export]]
Rcpp::List batch_estimator(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    double lam,
    arma::mat B_start,
    arma::mat Mu_start,
    arma::mat tau_start,
    int itermax,
    double stop_eps
    ){
    // initialization
    arma::rowvec Pi(K);
    Pi.fill(1.0 / K);
    arma::mat B(K,K), Mu(K,K), S(m,K);
    arma::mat tau(m,K);
    for (int k = 0; k < K; k++) {
        for (int l=0; l < K; l++) {
            B(k,l) = myrunif();
            Mu(k,l) = myrunif();
        }
    }
    // B.fill(0.5), Mu.fill(0.5);
    S.fill(1.0/K);
    tau.fill(1.0/K);
    //B = B_start, Mu = Mu_start;
    // for (int i = 0; i < m; i++) {
    //     arma::rowvec tt(K);
    //     for (int k = 0; k < K; k++) {
    //         tt(k) = myrunif();
    //     }
    //     tt = tt / sum(tt);
    //     tau.row(i) = tt;
    // }
    //tau = tau_start;

    int nall = alltimes.n_rows;
    int ncol = alltimes.n_cols;
    //int trunc_pos = 0, start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
    //int N = floor(T / dT);
    //int nsave = floor(5.0 / dT);
    //queue<double> trunc_pos_queue;
    
    //arma::vec elbo_vec(N);
    //double elbo = 0;
    //arma::mat prevdata;

    //double Tn, t_current, t_start, eta;
    //arma::rowvec event; 
    //arma::mat truncdata;
    Rcpp::List paralist;
    unordered_map<string, arma::vec> datamap;
    datamap = transfer(alltimes);


    double gap = 2147483647;
    double eta = 1.0/nall * (K * K);

    double t_start = 0.0, Tn;
    Tn = alltimes(nall - 1, ncol - 1);

    double trunc_length = 5.0;
    for (int iter = 0; iter < itermax; iter++) {
        eta = 1.0/nall * (K * K) /(iter + 1.0);
        // S.fill(0.0);
        paralist = update_lam_trunc(tau, Mu, B, Pi, S, datamap,
                                    t_start,
                                    Tn,
                                    m,
                                    K,
                                    A,
                                    lam,
                                    eta,
                                    trunc_length);
        arma::mat tau_new = paralist["tau"];
        arma::mat Mu_new = paralist["Mu"];
        arma::mat B_new = paralist["B"];
        arma::mat S_new = paralist["S"];
        arma::rowvec Pi_new = paralist["Pi"];
        double lam_new = paralist["lam"];
        gap = max(abs(Mu - Mu_new).max(), abs(B - B_new).max());
        tau = tau_new; 
        Mu = Mu_new, B = B_new, S = S_new, Pi = Pi_new;
        lam = lam_new;
        Rprintf("iter: %d \n", iter); 
        // B.print();
        // Mu.print();
        // Pi.print();
        Rprintf("lam: %2.3f", lam);
        Rprintf("gap: %2.3f", gap);
        Rprintf("=============\n");
        if (gap < stop_eps){
            break;
        }
    }

    return Rcpp::List::create(
                          Rcpp::Named("Mu") = Mu,
                          Rcpp::Named("B") = B,
                          Rcpp::Named("Pi") = Pi,
                          Rcpp::Named("lam") = lam,
                          Rcpp::Named("tau") = tau
                          );
}


// log-likelihood
double loglik_Hak(
    unordered_map<string, arma::vec> datamap,
    double t_start,
    double t_end,
    arma::vec Z,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    Rcpp::List A,
    double lam,
    int m,
    int K
    ){
    double loglik = 0;

    int l, k, n_edge;
    arma::vec timevec;
    string key;
    int ln,n;
    double intensity;
    unordered_map<string, arma::vec>:: iterator itr; 

    // arma::mat part1(K,K), part2(K,K), Lambda(K,K);
    double part1, part2, Lambda;

    double trunc_length = 5.0;
    double R = trunc_length / lam;

    for (itr = datamap.begin(); itr != datamap.end(); itr++) {
        key = itr->first;
        timevec = itr->second;
        arma::vec index = split(key);
        int i = (int) index(0), j = (int) index(1);
        ln = timevec.n_elem;
        k = Z(i), l = Z(j);

        part1 = 0.0, part2 = 0.0;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = 0.0;
                Lambda = 0.0;
                for (int n1 = n-1; n1 >= 0; n1--) {
                    double t1 = timevec(n1);
                    if (t_current - t1 > R)
                        break;
                    intensity += trigger(t1, t_current, lam);
                }

                Lambda = Mu(k,l) + B(k,l) * intensity;
                part1 = part1 + log(Lambda);
                part2 = part2 + B(k,l) * integral(t_current, t_end, lam);
            } else {
                part2 = part2 + B(k,l) * integral2(t_current, t_start, t_end, lam);
            }
        }

        loglik += (part1 - part2);
    }

    // mu * T
    for (int i = 0; i < m; i++) {
        arma::rowvec edge = A[i];
        n_edge = edge.n_elem;
        for (int p = 0; p < n_edge; p++) {
            int j = (int) edge(p);
            k = Z(i), l = Z(j);
            loglik -= Mu(k,l)* (t_end - t_start);
        }
        
    }

    // Pi
    for (int i = 0; i < m; i++) {
        k = Z(i);
        loglik += log(Pi(k) + eps);
    }

    return loglik;
}


// [[Rcpp::export]]
double get_loglik_Hak(
    arma::mat alltimes,
    double t_start,
    double t_end,
    arma::vec Z,
    arma::mat Mu,
    arma::mat B,
    arma::rowvec Pi,
    Rcpp::List A,
    double lam,
    int m,
    int K
    ){
    unordered_map<string, arma::vec> datamap = transfer(alltimes);
    double loglik = loglik_Hak(datamap, t_start, t_end, Z, Mu, B, Pi, A, lam, m, K);
    return loglik;
}
