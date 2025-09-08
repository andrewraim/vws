#ifndef RCPP_RECT_H
#define RCPP_RECT_H

#include <Rcpp.h>

//' Sequence
//'
//' @export
// [[Rcpp::export(name = "seq")]]
Rcpp::NumericVector rcpp_seq(double lo, double hi, unsigned int N, bool endpoints = false);

#endif

