#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <math.h>
#include <RcppArmadillo.h>
#include <vector>
#include <stdio.h>
#include <random>
#include <queue>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// ===============================
// miscellaneous functions
// ===============================

// append the second vecor to the first vector
arma::vec vecadd(
	arma::vec x,
	arma::vec y
	){
	int n1 = x.n_elem, n2 = y.n_elem;
	if (n2 == 0) { return x;}
	if (n1 == 0) { return y;}
	x.resize(n1 + n2);
	x.tail(n2) = y;
	return x;
}

// calculate number of events
int cal_event(unordered_map<string, arma::vec> datamap){
	int num = 0;
	unordered_map<string, arma::vec>:: iterator itr; 
	arma::vec event;
    for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
         event = itr->second;
         num += event.n_elem; 
     } 
     return num;
}

// calculate number of edges
// [[Rcpp::export]]
int cal_edge(Rcpp::List A){
	int m = A.size(), count = 0;
	for (int i = 0; i < m; i++) {
		arma::rowvec edge = A[i];
		count += edge.n_elem;
	}
	return count;
}


// 


// =============================
// Transfer matrix form into hash map form
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



// split string to two numbers separated by ","
// [[Rcpp::export]]
arma::vec split(string x){
	int pos = x.find(",");
	int l = x.length();
	arma::vec index(2);
	index(0) = stoi(x.substr(0,pos));
	index(1) = stoi(x.substr(pos + 1, l - pos - 1));
	return index;
}





// =====================================
// Update function for Hawkes
// =====================================



#endif
