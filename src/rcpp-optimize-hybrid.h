#ifndef RCPP_OPTIMIZE_HYBRID_H
#define RCPP_OPTIMIZE_HYBRID_H

#include <Rcpp.h>

//' Optimize Hybrid
//'
//' @param f TBD
//' @param init TBD
//' @param lower TBD
//' @param upper TBD
//' @param maximize TBD
//' @param maxiter TBD
//'
//' @export
// [[Rcpp::export(name = "optimize_hybrid")]]
Rcpp::List optimize_hybrid_rcpp(const Rcpp::Function& f, double init, double lower,
	double upper, bool maximize, unsigned maxiter = 100000);

#endif

