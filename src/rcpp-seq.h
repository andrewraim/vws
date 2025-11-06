#ifndef RCPP_SEQ_H
#define RCPP_SEQ_H

#include <Rcpp.h>

//' Produce a sequence of knots
//'
//' Produce knots which define \eqn{N} equally-spaced intervals between
//' (finite) endpoints `lo` and `hi`.
//'
//' @param lo Left endpoint; must be finite.
//' @param hi Right endpoint; must be finite.
//' @param N Number of desired intervals.
//' @param endpoints logical; if `TRUE`, include the endpoints.
//'
//' @returns
//' A vector that represents a sequence of knots. If `endpoints = TRUE`, it
//' contains \eqn{N+1} evenly-spaced knots that represent \eqn{N} regions with
//' endpoints included. If `endpoints = FALSE`, the endpoints are excluded.
//'
//' @examples
//' seq_knots(0, 1, N = 5)
//' seq_knots(0, 1, N = 5, endpoints = TRUE)
//'
//' # Trivial case: make endpoints for just one interval
//' seq_knots(0, 1, N = 1)
//' seq_knots(0, 1, N = 1, endpoints = TRUE)
//'
//' # The following calls throw errors
//' tryCatch({
//'   seq_knots(0, 1, N = 0)
//' }, error = function(e) { print(e) })
//' tryCatch({
//'   seq_knots(0, Inf, N = 5)
//' }, error = function(e) { print(e) })
//' tryCatch({
//'   seq_knots(-Inf, 1, N = 5)
//' }, error = function(e) { print(e) })
//'
//' @export
// [[Rcpp::export(name = "seq_knots")]]
 Rcpp::NumericVector rcpp_seq(double lo, double hi, unsigned int N,
	bool endpoints = false);

#endif
