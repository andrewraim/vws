#ifndef VWS_UNIVARIATE_HELPER_H
#define VWS_UNIVARIATE_HELPER_H

#include <Rcpp.h>

namespace vws {

//' Univariate Distribution Helper
//'
//' A container to support VWS with univariate distributions.
//'
//' @param d Density function.
//' @param p Cumulative distribution function.
//' @param q Quantile function.
//' @param s Indicator function that returns \code{TRUE} when the
//' argument is in the support of the distribution; otherwise returns
//' \code{FALSE}.
//'
//' @details
//' The specified functions must support the following interfaces.
//' \describe{
//' \item{\code{d(x, log = FALSE)}}{}
//' \item{\code{p(q, lower.tail = TRUE, log.p = FALSE)}}{}
//' \item{\code{q(p, lower.tail = TRUE, log.p = FALSE)}}{}
//' \item{\code{s(x)}}{}
//' }
//'
//' Arguments to these functions are interpreted as usual for the \code{stats}
//' package. (Other arguments will be ignored if present).
//' \describe{
//' \item{\code{x}, \code{q}}{vector of quantiles.}{}
//' \item{\code{log}, \code{log.p}}{logical; if \code{TRUE}, probabilities
//' \code{p} are given as \eqn{\log(p)}.}{}
//' \item{\code{lower.tail}}{logical; if \code{TRUE}, probabilities are
//' \eqn{\text{P}(X \leq x)}; otherwise \eqn{\text{P}(X > x)}}{}
//' }
//'
//' @examples
//' helper = univariate_helper(
//'     d = function(x, log = FALSE) dpois(x, 10, log),
//'     p = function(q, lower.tail = TRUE, log.p = FALSE) {
//'         ppois(q, 10, lower.tail, log.p)
//'     },
//'     q = function(p, lower.tail = TRUE, log.p = FALSE) {
//'         qpois(p, 10, lower.tail, log.p)
//'     },
//'     s = function(x) { is.numeric(x) }
//' )
//'
//' helper$d(5)
//' helper$p(10)
//' helper$q(0.025)
//'
//' @name univariate_helper
//' @export
template <class T>
class UnivariateHelper
{
public:
	virtual double pdf(T x, bool log = false) const = 0;
	virtual double cdf(T q, bool lower = true, bool log = false) const = 0;
	virtual double quantile(T p, bool lower = true, bool log = false) const = 0;
	virtual bool supp(T x) const = 0;
};

}

#endif
