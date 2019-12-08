#include "onlineblock.h"

unordered_map<string, arma::mat> w_normaliztion(arma::mat W1, arma::mat W2){
	int K = W1.n_cols;
	double r1, r2, r0;
	unordered_map<string, arma::mat> res;
	for(int k = 0; k < K; k++) {
		r1 = sum(W1.col(k));
		r2 = sum(W2.col(k));
		r0 = sqrt(r1 * r2);
		W1.col(k) = W1.col(k) / r1 * r0;
		W2.col(k) = W2.col(k) / r2 * r0;
	}
	res["W1"] = W1;
	res["W2"] = W2;
	return res;
}


Rcpp::List update_ccrm(
    arma::mat W1,
    arma::mat W2,
    double b,
    unordered_map<string, std::deque<double>> datamap,
    double t_start,
    double t_end,
    int m,
    int K,
    Rcpp::List A,
    double lam,
    double eta
    ){
    arma::mat P1_w1(m,K), P2_w1(m,K), P1_w2(m,K), P2_w2(m,K);
    double  P1_b, P2_b;
    arma::mat P1_w1_tp(m,K), P2_w1_tp(m,K), P1_w2_tp(m,K), P2_w2_tp(m,K);
    double P1_b_tp, P2_b_tp, Lambda = 0.0;
    P1_w1.fill(0.0), P2_w1.fill(0.0), P1_w2.fill(0.0), P2_w2.fill(0.0);
    P1_b = P2_b = 0.0;
    double lam_store;
    double mu;

    int k;
    //int l, k;
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
        

        mu = 0.0;
       	for (k = 0; k < K; k++) {
            P2_w1(i,k) = P2_w1(i,k) + W2(j,k) * (t_end - t_start);
            P2_w2(j,k) = P2_w2(j,k) + W1(i,k) * (t_end - t_start);
            mu += W1(i,k) * W2(j,k);
        }

        P1_w1_tp.fill(0.0), P2_w1_tp.fill(0.0), P1_w2_tp.fill(0.0), P2_w2_tp.fill(0.0);
        P1_b_tp = 0.0, P2_b_tp = 0.0;
        lam_store = 0.0;
        //P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda = eps; // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_b_tp = P2_b_tp + integral(t_current, t_end, lam);
                Lambda += mu + b * intensity;

                lam_store = lam_store + b * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - b * ((t_end - t_current) * exp(-lam*(t_end - t_current)));

                P1_b_tp += intensity / Lambda;
                for (k = 0; k < K; k++) {
                	P1_w1_tp(i,k) += W2(j,k) / Lambda;
                	P1_w2_tp(j,k) += W1(i,k) / Lambda;
                }
            } else {
                P2_b_tp = P2_b_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + b * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        grad_lam += lam_store;

        /*
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
            }
        }
        */

        for (k = 0; k < K; k++) {
        	P1_w1(i,k) += P1_w1_tp(i,k);
        	P1_w2(j,k) += P1_w2_tp(j,k);
        }
        P1_b += P1_b_tp;
        P2_b += P2_b_tp;

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
    double grad_b = P1_b - P2_b;
    //grad_B.print();
    arma::mat grad_w1 = P1_w1 - P2_w1;
    arma::mat grad_w2 = P1_w2 - P2_w2;
    //grad_mu.print();
    double b_new = b + eta * grad_b;
    //printf("B new is: \n");   
    //B_new.print();
    arma::mat W1_new = W1 + eta * m * grad_w1;
    arma::mat W2_new = W2 + eta * m * grad_w2;
    //printf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    if (b_new <= 0.0) 
    	b_new = b / 2.0;
    else if (b_new >= 2*b)
    	b_new = 2.0 * b;

    for (int i = 0; i < m; i++) {
    	for(k = 0; k < K; k++) {
    		if (W1_new(i,k) <= 0.0) 
                W1_new(i,k) = W1(i,k) / 2.0;
            else if (W1_new(i,k) > 2 * W1(i,k))
                W1_new(i,k) = W1(i,k) * 2.0;

            if (W2_new(i,k) <= 0.0) 
                W2_new(i,k) = W2(i,k) / 2.0;
            else if (W2_new(i,k) > 2 * W2(i,k))
                W2_new(i,k) = W2(i,k) * 2.0;
    	}
    }

    
    // normaliztion
    unordered_map<string, arma::mat> res = w_normaliztion(W1_new, W2_new);
    W1_new = res["W1"];
    W2_new = res["W2"];
	

    double lam_new = lam + eta * grad_lam;
    if (lam_new > 5*lam) {
        lam_new = 5 * lam;
    } else if (lam_new <= 0.0) {
        lam_new = lam/2.0;
    }


    return Rcpp::List::create(Rcpp::Named("W1") = W1_new,
                          Rcpp::Named("W2") = W2_new,
                          Rcpp::Named("b") = b_new,
                          Rcpp::Named("lam") = lam_new);
}



// [[Rcpp::export]]
Rcpp::List ccrm_estimator(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    double lam,
    arma::mat W1_start,
    arma::mat W2_start,
    double b_start,
    bool is_elbo = false
    ){

    // create empty map
    unordered_map<string, std::deque<double>> datamap;
    datamap = transfer_create(A, m);

    // initialization
    arma::mat W1(m,K), W2(m,K);
    double b = myrunif();
    for (int i = 0; i < m; i++) {
        for (int k = 0; k < K; k++) {
            W1(i,k) = 2*myrunif();
            W2(i,k) = 2*myrunif();
        }
    }
    // W1.fill(1.0);
    //b = b_start, W1 = W1_start, W2 = W2_start;

    int nall = alltimes.n_rows;
    int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
    int N = floor(T / dT);

    double R = 5.0;
    
    arma::vec elbo_vec(N);
    //double elbo = 0.0;
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
        // datamap = transfer(truncdata);
        datamap = transfer_eff(datamap, truncdata, R);

        t_start = Tn - dT;
        ln_curr = end_pos;
        n_t = ln_curr - ln_prev;
        // eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
        eta = 1.0/sqrt(1 + n/10.0)/n_t;
        paralist = update_ccrm(W1, W2, b, datamap, t_start, Tn, m, K, A, lam, eta);
        // paralist = update_lam_eff(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
        arma::mat W1_new = paralist["W1"], W2_new = paralist["W2"];
        double lam_new = paralist["lam"], b_new = paralist["b"]; 
        W1 = W1_new;
        W2 = W2_new;
        lam = lam_new, b = b_new;
        start_pos = curr_pos;
        ln_prev = ln_curr;
        printf("iter: %d; number: %d \n", n, n_t); 
        printf("b: %2.3f", b);
        printf("lam: %2.3f", lam);

        if (is_elbo) {
        	//elbo = 0.0;
            //prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
            //elbo = get_elbo_Hak(prevdata, 0.0, Tn, tau, Mu, B, Pi, A, lam, m, K);
            //elbo_vec(n) = elbo / ln_curr;
        }

        //S.print();
        printf("=============\n");
    }

    return Rcpp::List::create(
                          Rcpp::Named("W1") = W1,
                          Rcpp::Named("W2") = W2,
                          Rcpp::Named("b") = b,
                          Rcpp::Named("lam") = lam,
                          Rcpp::Named("elbo") = elbo_vec);
}



// --- test ----

Rcpp::List update_test(
    arma::mat W,
    double b,
    unordered_map<string, std::deque<double>> datamap,
    double t_start,
    double t_end,
    int m,
    int K,
    Rcpp::List A,
    double lam,
    double eta
    ){
    arma::mat P1_w(m,m), P2_w(m,m);
    double  P1_b, P2_b;
    arma::mat P1_w_tp(m,m), P2_w_tp(m,m);
    double P1_b_tp, P2_b_tp, Lambda = 0.0;
    P1_w.fill(0.0), P2_w.fill(0.0);
    P1_b = P2_b = 0.0;
    double lam_store;
    double mu;

    //int l, k;
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
        

        mu = 0.0;
       	mu = W(i,j);
       	P2_w(i,j) = P2_w(i,j) + t_end - t_start;

        P1_w_tp.fill(0.0), P2_w_tp.fill(0.0);
        P1_b_tp = 0.0, P2_b_tp = 0.0;
        lam_store = 0.0;
        //P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda = eps; // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_b_tp = P2_b_tp + integral(t_current, t_end, lam);
                Lambda += mu + b * intensity;

                lam_store = lam_store + b * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - b * ((t_end - t_current) * exp(-lam*(t_end - t_current)));

                P1_b_tp += intensity / Lambda;
                P1_w(i,j) += 1/Lambda;
            } else {
                P2_b_tp = P2_b_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + b * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        grad_lam += lam_store;

        /*
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
            }
        }
        */


        P1_w(i,j) += P1_w_tp(i,j);
        P1_b += P1_b_tp;
        P2_b += P2_b_tp;
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
    double grad_b = P1_b - P2_b;
    //grad_B.print();
    arma::mat grad_w = P1_w - P2_w;
    //grad_mu.print();
    double b_new = b + eta * grad_b;
    //printf("B new is: \n");   
    //B_new.print();
    arma::mat W_new = W + eta * m * grad_w;
    //printf("Mu new is: \n");
    //Mu_new.print();

    // handle negative values and large gradient
    if (b_new <= 0.0) 
    	b_new = b / 2.0;
    else if (b_new >= 2*b)
    	b_new = 2.0 * b;

    for (int i = 0; i < m; i++) {
    	for (int j = 0; j < m; j++){
    		if (i != j){
    			if (W_new(i,j) <= 0.0)
    				W_new(i,j) = W(i,j)/2.0;
    			else if (W_new(i,j) >= 2*W(i,j))
    				W_new(i,j) = W(i,j) * 2.0;
    		}
    	}
    }
    

    double lam_new = lam + eta * grad_lam;
    if (lam_new > 5*lam) {
        lam_new = 5 * lam;
    } else if (lam_new <= 0.0) {
        lam_new = lam/2.0;
    }


    return Rcpp::List::create(Rcpp::Named("W") = W_new,
                          Rcpp::Named("b") = b_new,
                          Rcpp::Named("lam") = lam_new);
}



// [[Rcpp::export]]
Rcpp::List test_estimator(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    double lam,
    arma::mat W_start,
    double b_start,
    bool is_elbo = false
    ){

    // create empty map
    unordered_map<string, std::deque<double>> datamap;
    datamap = transfer_create(A, m);

    // initialization
    arma::mat W(m,m);
    W.fill(0.0);
    double b = myrunif();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
        	if (i == j)
        		continue;
            W(i,j) = myrunif();
        }
    }
    //b = b_start, W1 = W1_start, W2 = W2_start;

    int nall = alltimes.n_rows;
    int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
    int N = floor(T / dT);

    double R = 5.0;
    
    arma::vec elbo_vec(N);
    //double elbo = 0.0;
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
        // datamap = transfer(truncdata);
        datamap = transfer_eff(datamap, truncdata, R);

        t_start = Tn - dT;
        ln_curr = end_pos;
        n_t = ln_curr - ln_prev;
        // eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
        eta = 1.0/sqrt(1 + n/10.0)/n_t;
        //printf("here 1 \n");
        paralist = update_test(W, b, datamap, t_start, Tn, m, K, A, lam, eta);
        // paralist = update_lam_eff(tau, Mu, B, Pi, S, datamap, t_start, Tn, m, K, A, lam, eta);
        arma::mat W_new = paralist["W"];
        double lam_new = paralist["lam"], b_new = paralist["b"]; 
        W = W_new;
        lam = lam_new, b = b_new;
        start_pos = curr_pos;
        ln_prev = ln_curr;
        printf("iter: %d; number: %d \n", n, n_t); 
        printf("b: %2.3f", b);
        printf("lam: %2.3f", lam);

        if (is_elbo) {
        	//elbo = 0.0;
            //prevdata = alltimes.rows(0, end_pos - 1); // head_rows()
            //elbo = get_elbo_Hak(prevdata, 0.0, Tn, tau, Mu, B, Pi, A, lam, m, K);
            //elbo_vec(n) = elbo / ln_curr;
        }

        //S.print();
        printf("=============\n");
    }

    return Rcpp::List::create(
                          Rcpp::Named("W") = W,
                          Rcpp::Named("b") = b,
                          Rcpp::Named("lam") = lam,
                          Rcpp::Named("elbo") = elbo_vec);
}

Rcpp::List compute_grad_ccrm(
    arma::mat W1,
    arma::mat W2,
    double b,
    unordered_map<string, std::deque<double>> datamap,
    double t_start,
    double t_end,
    int m,
    int K,
    Rcpp::List A,
    double lam,
    double eta
    ){
    arma::mat P1_w1(m,K), P2_w1(m,K), P1_w2(m,K), P2_w2(m,K);
    double  P1_b, P2_b;
    arma::mat P1_w1_tp(m,K), P2_w1_tp(m,K), P1_w2_tp(m,K), P2_w2_tp(m,K);
    double P1_b_tp, P2_b_tp, Lambda = 0.0;
    P1_w1.fill(0.0), P2_w1.fill(0.0), P1_w2.fill(0.0), P2_w2.fill(0.0);
    P1_b = P2_b = 0.0;
    double lam_store;
    double mu;

    int k;
    //int l, k;
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
        

        mu = 0.0;
       	for (k = 0; k < K; k++) {
            P2_w1(i,k) = P2_w1(i,k) + W2(j,k) * (t_end - t_start);
            P2_w2(j,k) = P2_w2(j,k) + W1(i,k) * (t_end - t_start);
            mu += W1(i,k) * W2(j,k);
        }

        P1_w1_tp.fill(0.0), P2_w1_tp.fill(0.0), P1_w2_tp.fill(0.0), P2_w2_tp.fill(0.0);
        P1_b_tp = 0.0, P2_b_tp = 0.0;
        lam_store = 0.0;
        //P2_mu_tp = P2_mu_tp + t_end - t_start;
        ln = timevec.n_elem;
        for (n = ln - 1; n >= 0; n--){
            double t_current = timevec(n);
            if (t_current > t_start) {
                intensity = eps, intensity_lam1 = eps, intensity_lam2 = eps;
                Lambda = eps; // store temporary intensity values
                for (int n1 = 0; n1 < n; n1++) {
                    double t1 = timevec(n1);
                    intensity += trigger(t1, t_current, lam);
                    intensity_lam1 += trigger_lam(t1, t_current, lam);
                    intensity_lam2 += (t_current - t1) * trigger(t1, t_current, lam);
                }
                P2_b_tp = P2_b_tp + integral(t_current, t_end, lam);
                Lambda += mu + b * intensity;

                lam_store = lam_store + b * (intensity_lam1 - intensity_lam2) / Lambda;
                lam_store = lam_store - b * ((t_end - t_current) * exp(-lam*(t_end - t_current)));

                P1_b_tp += intensity / Lambda;
                for (k = 0; k < K; k++) {
                	P1_w1_tp(i,k) += W2(j,k) / Lambda;
                	P1_w2_tp(j,k) += W1(i,k) / Lambda;
                }
            } else {
                P2_b_tp = P2_b_tp + integral2(t_current, t_start, t_end, lam);
                lam_store = lam_store + b * ((t_start - t_current) * exp(-lam*(t_start - t_current)) - (t_end - t_current) * exp(-lam*(t_end - t_current)));
            }
        }

        grad_lam += lam_store;

        /*
        for (k = 0; k < K; k++) {
            for (l = 0; l < K; l++) {
                grad_lam += tau(i,k) * tau(j,l) * lam_store(k,l);
            }
        }
        */

        for (k = 0; k < K; k++) {
        	P1_w1(i,k) += P1_w1_tp(i,k);
        	P1_w2(j,k) += P1_w2_tp(j,k);
        	P1_b += P1_b_tp;
        	P2_b += P2_b_tp;
        }

    } 

    // update parameters
    double grad_b = P1_b - P2_b;
    //grad_B.print();
    arma::mat grad_w1 = P1_w1 - P2_w1;
    arma::mat grad_w2 = P1_w2 - P2_w2;
    //grad_mu.print();
    


    return Rcpp::List::create(Rcpp::Named("W1") = grad_w1,
                          Rcpp::Named("W2") = grad_w2,
                          Rcpp::Named("b") = grad_b,
                          Rcpp::Named("lam") = grad_lam);
}


// [[Rcpp::export]]
void get_grad_ccrm(
    arma::mat alltimes,
    Rcpp::List A,
    int m,
    int K,
    double T,
    double dT,
    double lam,
    arma::mat W1_start,
    arma::mat W2_start,
    double b_start
    ){

    // create empty map
    unordered_map<string, std::deque<double>> datamap;
    datamap = transfer_create(A, m);

    // initialization
    arma::mat W1(m,K), W2(m,K);
    double b;
    b = b_start, W1 = W1_start, W2 = W2_start;

    int nall = alltimes.n_rows;
    int start_pos = 0, curr_pos = 0, end_pos = 0, ln_prev = 0, ln_curr, n_t;
    int N = floor(T / dT);

    double R = 20.0;
    
    arma::vec elbo_vec(N);
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
        // datamap = transfer(truncdata);
        datamap = transfer_eff(datamap, truncdata, R);

        t_start = Tn - dT;
        ln_curr = end_pos;
        n_t = ln_curr - ln_prev;
        // eta = 1.0/sqrt(1 + n/10.0)/n_t * (K * K);
        eta = 1.0/sqrt(1 + n/10.0)/n_t;
        // paralist = update_ccrm(W1, W2, b, datamap, t_start, Tn, m, K, A, lam, eta);
        paralist = compute_grad_ccrm(W1, W2, b, datamap, t_start, Tn, m, K, A, lam, eta);
        arma::mat grad_w1 = paralist["W1"], grad_w2 = paralist["W2"];
        double grad_lam = paralist["lam"], grad_b = paralist["b"]; 

        start_pos = curr_pos;
        ln_prev = ln_curr;
        printf("iter: %d; number: %d \n", n, n_t); 
        grad_w1.print();
        grad_w2.print();
        printf("b: %2.3f \n", grad_b);
        printf("lam: %2.3f \n", grad_lam);

        //S.print();
        printf("=============\n");
    }

}
