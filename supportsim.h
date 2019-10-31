#include "onlineblock.h"


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