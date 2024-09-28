#ifndef RCPP_GUMBEL_H
#define RCPP_GUMBEL_H

#include <Rcpp.h>

//' Gumbel Distribution
//'
//' Functions for the Gumbel distribution
//'
//' @param n Number of draws
//' @param x Vector of quantiles
//' @param p Vector of probabilities
//' @param q Vector of quantiles
//' @param mu Location parameter
//' @param sigma Scale parameter
//' @param lower Logical; if `TRUE` (default), probabilities are
//' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
//' @param log Logical; if `TRUE`, probabilities p are given as \eqn{log(p)}
//'
//' @return A vector of draws
//'
//' @name Gumbel
//' @export
// [[Rcpp::export(name = "r_gumbel")]]
 Rcpp::NumericVector r_gumbel_rcpp(unsigned int n, double mu = 0, double sigma = 1);

//' @name Gumbel
//' @export
// [[Rcpp::export(name = "d_gumbel")]]
Rcpp::NumericVector d_gumbel_rcpp(const Rcpp::NumericVector& x, double mu = 0,
	double sigma = 1, bool log = false);

//' @name Gumbel
//' @export
// [[Rcpp::export(name = "p_gumbel")]]
Rcpp::NumericVector p_gumbel_rcpp(const Rcpp::NumericVector& q, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false);

//' @name Gumbel
//' @export
// [[Rcpp::export(name = "q_gumbel")]]
Rcpp::NumericVector q_gumbel_rcpp(const Rcpp::NumericVector& p, double mu = 0,
	double sigma = 1, bool lower = true, bool log = false);

#endif
