#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <math.h>
#include <RcppArmadillo.h>
#include <vector>
#include <stdio.h>
#include <random>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



// =============================
// Transfer matrix forn into hash map form
// =============================

// [[Rcpp::export]]
unordered_map<string, arma::vec> transfer(arma::mat alltimes){
	unordered_map<string, arma::vec> datamap;
	int N = alltimes.n_rows;
	arma::rowvec event;
	int i,j,ln;
	double time;
	string key;
	
	for (int n = 0; n < N; n++) {
		event = alltimes.row(n);
		i = event[0], j = event[1], time = event[2];
		key = to_string(i) + "," + to_string(j);
		if (datamap.find(key) == datamap.end()) {
			// key not found 
			arma::vec timevec(1);
			timevec(0) = time;
			datamap[key] = timevec;
		} else {
			// key found
			arma::vec timevec = datamap[key];
			ln = timevec.n_elem;
			timevec.resize(ln + 1);
			timevec(ln) = time;
			datamap[key] = timevec;
		}
	}
	return datamap;
}

// this function cannot be exported, since Rcpp cannot convert unordered_map
void mapprint(unordered_map<string, arma::vec> datamap){
	unordered_map<string, arma::vec>:: iterator itr; 
    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
        cout << itr->first << "  " << itr->second << endl; 
     } 
}

// [[Rcpp::export]]
void testmap(arma::mat alltimes){
	// print map for check
	unordered_map<string, arma::vec> datamap = transfer(alltimes);
	mapprint(datamap);
}

// =====================================
// Update function for homo Poisson
// =====================================


// =====================================
// Update function for Hawkes
// =====================================

#endif
