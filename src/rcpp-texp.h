#ifndef RCPP_TEXP_H
#define RCPP_TEXP_H

#include <Rcpp.h>

//' Truncated Exponential
//'
//' @name TExp
//' @export
// [[Rcpp::export(name = "n_texp")]]
double n_texp_rcpp(double kappa, double lo, double hi, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "integrate_texp")]]
double integrate_texp_rcpp(double a, double b, double kappa, double lo,
	double hi, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "d_texp")]]
Rcpp::NumericVector d_texp_rcpp(const Rcpp::NumericVector& x, double kappa,
	double lo, double hi, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "p_texp")]]
Rcpp::NumericVector p_texp_rcpp(const Rcpp::NumericVector& q, double kappa,
	double lo, double hi, bool lower = true, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "q_texp")]]
Rcpp::NumericVector q_texp_rcpp(const Rcpp::NumericVector& p, double kappa,
	double lo, double hi, bool lower = true, bool log = false);

//' @name TExp
//' @export
// [[Rcpp::export(name = "r_texp")]]
Rcpp::NumericVector r_texp_rcpp(unsigned int n, double kappa, double lo,
	double hi);

//' @name TExp
//' @export
// [[Rcpp::export(name = "mgf_texp")]]
Rcpp::NumericVector mgf_texp_rcpp(const Rcpp::NumericVector& s, double kappa,
	double lo, double hi, bool log = false);

#endif
