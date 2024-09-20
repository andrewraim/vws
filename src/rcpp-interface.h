#ifndef RCPP_INTERFACE_H
#define RCPP_INTERFACE_H

#include <Rcpp.h>
#include "vws.h"

//' Rectangular transformation
//'
//' A transformation from Euclidean to a rectangle in \eqn{\mathbb{R}^n}, and
//' its inverse transformation.
//'
//' @param x A point in \eqn{\mathbb{R}^{d}}.
//' @param z A point in the rectangle \eqn{[a_1,b_1] \times \cdots \times [a_d,b_d]}.
//' @param a A vector \eqn{(a_1, \ldots, a_d)}, Elements may be \code{-Inf}.
//' @param b A vector \eqn{(b_1, \ldots, b_d)}, Elements may be \code{Inf}.
//'
//' @examples
//' x = seq(-1, 1, length.out = 3)
//' a = rep(0, 3)
//' b = rep(1, 3)
//' z = inv_rect(x, a, b)
//' rect(z, a, b)
//'
//' a = c(-Inf, 0, -Inf)
//' b = c(Inf, 1, Inf)
//' z = inv_rect(x, a, b)
//' rect(z, a, b)
//'
//' @name rect
//' @export
Rcpp::NumericVector rect_rcpp(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b);

//' @name rect
//' @export
Rcpp::NumericVector inv_rect_rcpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b);

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
