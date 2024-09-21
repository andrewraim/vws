#ifndef RCPP_INVGAMMA_H
#define RCPP_INVGAMMA_H

#include <Rcpp.h>

//' Inverse Gamma distribution
//'
//' @param n Number of draws.
//' @param x Vector of quantiles.
//' @param q Vector of quantiles.
//' @param p Vector of probabilities.
//' @param a Shape parameter.
//' @param b Rate parameter.
//' @param lower logical; if TRUE (default), probabilities are
//' \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
//' @param log If `TRUE`, return densities and probabilities on the log-scale.
//'
//' @return
//' `d_invgamma` computes the density, `r_invgamma` generates random deviates,
//' `p_invgamma` computes the CDF, and `q_invgamma` computes quantiles.
//'
//' @details
//' Note that `Rcpp::*gamma` and `R::*gamma` functions are both parameterized
//' by a scale parameter, which is the inverse of the rate.
//'
//' @name InverseGamma
// [[Rcpp::export(name = "r_invgamma")]]
Rcpp::NumericVector r_invgamma_rcpp(unsigned int n, double a, double b);

//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "d_invgamma")]]
Rcpp::NumericVector d_invgamma_rcpp(const Rcpp::NumericVector& x, double a,
	double b, bool log = false);

//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "p_invgamma")]]
Rcpp::NumericVector p_invgamma_rcpp(const Rcpp::NumericVector& q, double a,
	double b, bool lower = true, bool log = false);

//' @name InverseGamma
//' @export
// [[Rcpp::export(name = "q_invgamma")]]
Rcpp::NumericVector q_invgamma_rcpp(const Rcpp::NumericVector& p, double a,
	double b, bool lower = true, bool log = false);

#endif
