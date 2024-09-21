#include "RUnivariateHelper.h"

RUnivariateHelper::RUnivariateHelper(const Rcpp::Function& pdf,
	const Rcpp::Function& cdf, const Rcpp::Function& quantile, const Rcpp::Function& supp)
: UnivariateHelper(), _pdf(pdf), _cdf(cdf), _quantile(quantile), _supp(supp)
{
}

double RUnivariateHelper::pdf(double x, bool log) const {
	const Rcpp::NumericVector& out = _pdf(x, log);
	return out(0);
}

double RUnivariateHelper::cdf(double x, bool lower, bool log) const {
	const Rcpp::NumericVector& out = _cdf(x, lower, log);
	return out(0);
}

double RUnivariateHelper::quantile(double x, bool lower, bool log) const {
	const Rcpp::NumericVector& out = _quantile(x, lower, log);
	return out(0);
}

bool RUnivariateHelper::supp(double x) const {
	const Rcpp::LogicalVector& out = _supp(x);
	return out(0);
}
