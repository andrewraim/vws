#ifndef RCPP_SEQ_H
#define RCPP_SEQ_H

#include <Rcpp.h>

//' Sequence
//'
//' @export
Rcpp::NumericVector rcpp_seq(double lo, double hi, unsigned int N, bool endpoints = false);

#endif

