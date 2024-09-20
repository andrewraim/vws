#include "RUnivariateHelper.h"

RUnivariateHelper::RUnivariateHelper()
: _d(NULL), _p(NULL), _q(NULL), _s(NULL)
{
	Rprintf("In empty constructor of RUnivariateHelper\n");
}


RUnivariateHelper::RUnivariateHelper(const Rcpp::Function& d, const Rcpp::Function& p,
	const Rcpp::Function& q, const Rcpp::Function& s)
: _d(d), _p(p), _q(q), _s(s)
{
	Rprintf("In constructor of RUnivariateHelper\n");
}

double RUnivariateHelper::d(double x, bool log) const {
	Rprintf("RUnivariateHelper: called d\n");
	return 0;
	// const Rcpp::NumericVector& out = _d(x, log);
	// return out(0);
}

double RUnivariateHelper::p(double q, bool lower, bool log) const {
	Rprintf("RUnivariateHelper: called p\n");
	return 0;
	// const Rcpp::NumericVector& out = _p(q, lower, log);
	// return out(0);
}

double RUnivariateHelper::q(double p, bool lower, bool log) const {
	Rprintf("RUnivariateHelper: called q\n");
	return 0;
	// const Rcpp::NumericVector& out = _q(p, lower, log);
	// return out(0);
}

bool RUnivariateHelper::s(double x) const {
	Rprintf("RUnivariateHelper: called s\n");
	return 0;
	// const Rcpp::LogicalVector& out = _s(x);
	// return out(0);
}
