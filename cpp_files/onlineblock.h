#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <math.h>
#include <RcppArmadillo.h>
#include <vector>
#include <stdio.h>
#include <random>
#include <queue>
#include <deque>

#include <chrono>
using namespace std::chrono;

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// define some constants

const double eps = 0.000001;

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



// [[Rcpp::export]]
unordered_map<string, arma::vec> transfer2(arma::mat alltimes, Rcpp::List A, int m){
	unordered_map<string, arma::vec> datamap;
	int N = alltimes.n_rows;
	arma::rowvec event;
	int i,j,ln;
	double time;
	string key;
	int n_edge;
	
	for (int i = 0; i < m; i++) {
        arma::rowvec edge = A[i];
        n_edge = edge.n_elem;
        for (int p = 0; p < n_edge; p++){
            int j = (int) edge(p);
            key = to_string(i) + "," + to_string(j);
            arma::vec timevec;
            datamap[key] = timevec;
        }
    }

    arma::vec timetemp;
	for (int n = 0; n < N; n++) {
		event = alltimes.row(n);
		i = event[0], j = event[1], time = event[2];
		key = to_string(i) + "," + to_string(j);
		timetemp = datamap[key];
		ln = timetemp.n_elem;
		timetemp.resize(ln + 1);
		timetemp(ln) = time;
		datamap[key] = timetemp;
	}
	return datamap;
}



// create a empty map
// [[Rcpp::export]]
unordered_map<string, std::deque<double> > transfer_create(Rcpp::List A, int m){
	unordered_map<string, std::deque<double>> datamap;
	int i,j;
	string key;
	int n_edge;
	
	for (i = 0; i < m; i++) {
        arma::rowvec edge = A[i];
        n_edge = edge.n_elem;
        for (int p = 0; p < n_edge; p++){
            j = (int) edge(p);
            key = to_string(i) + "," + to_string(j);
            std::deque<double> timevec;
            datamap[key] = timevec;
        }
    }
	return datamap;
}


// trim queue, old data exceeding R will be deleted
std::deque<double> trim_queue(std::deque<double> data, double t_current, double R){
	data.push_back(t_current);
	if (data.size() == 1){
		return data;
	}
	double time = data.front();
	while (t_current - time > R) {
		data.pop_front();
		time = data.front();
	}
	return data;
}

// efficient way to trasnfer data to reduce the head cost

unordered_map<string, std::deque<double>> transfer_eff2(unordered_map<string, std::deque<double>> datamap, arma::mat newtimes, double R){
	int N = newtimes.n_rows;
	arma::rowvec event;
	int i,j;
	double time;
	string key;

    std::deque<double> timetemp;
	for (int n = 0; n < N; n++) {
		event = newtimes.row(n);
		i = event[0], j = event[1], time = event[2];
		key = to_string(i) + "," + to_string(j);
		timetemp = datamap[key];
		timetemp = trim_queue(timetemp, time, R);
		datamap[key] = timetemp;
	}
	return datamap;
}


// use reference 
void transfer_eff(unordered_map<string, std::deque<double>> &datamap, arma::mat newtimes, double R){
	int N = newtimes.n_rows;
	arma::rowvec event;
	int i,j;
	double time;
	string key;

    std::deque<double> timetemp;
	for (int n = 0; n < N; n++) {
		event = newtimes.row(n);
		i = event[0], j = event[1], time = event[2];
		key = to_string(i) + "," + to_string(j);
		timetemp = datamap[key];
		timetemp = trim_queue(timetemp, time, R);
		datamap[key] = timetemp;
	}
}



unordered_map<string, std::deque<double> > transfer_create_empty(){
	unordered_map<string, std::deque<double>> datamap;
	return datamap;
}

// do not keep all datamap
void transfer_dynamic(unordered_map<string, std::deque<double>> &datamap, arma::mat newtimes, double R, double t_current){
	int N = newtimes.n_rows;
	arma::rowvec event;
	int i,j;
	double time;
	string key;

    std::deque<double> timetemp;
    std::unordered_map<std::string, std::deque<double>>::iterator got;
	for (int n = 0; n < N; n++) {
		event = newtimes.row(n);
		i = event[0], j = event[1], time = event[2];
		key = to_string(i) + "," + to_string(j);

		got = datamap.find(key);
		if (got == datamap.end()){
			// if key in data map
			timetemp.clear();
			timetemp.push_back(time);
			datamap[key] = timetemp;
		} else {
			// else if key is not in data map
			timetemp = datamap[key];
			timetemp = trim_queue(timetemp, time, R);
			datamap[key] = timetemp;			
		}				
	}
	// loop over datamap to check whether need to throw old pairs
	unordered_map<string, std::deque<double>>:: iterator itr; 
	std::deque<double> timeque;
	for (itr = datamap.begin(); itr != datamap.end(); itr++) 
    { 
        // type itr->first stores the key part  and 
        // itr->second stroes the value part 
        key = itr->first;
        timeque = itr->second;
        if (timeque.back() < t_current - R) {
        	// remove 
        	//cout<<"In here?"<<endl;
        	cout<<key<<endl;
        	datamap.erase(key);
        }
    }
}



// convert back to arma::vec
arma::vec convert_deque(std::deque<double> mydeque){
	int ln = mydeque.size();
	int count = 0;
	arma::vec myvec(ln);
	std::deque<double>::iterator it = mydeque.begin();

  	while (it != mydeque.end()){
  		myvec(count) = *it;
    	it++;
    	count++; 
  	}
  	return myvec;
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
// Update function for homo Poisson
// =====================================


// =====================================
// Update function for Hawkes
// =====================================


#endif
