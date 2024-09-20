#ifndef R_UNIVARIATE_HELPER_H
#define R_UNIVARIATE_HELPER_H

#include "vws.h"

class RUnivariateHelper : public vws::UnivariateHelper<double>
{
public:
	RUnivariateHelper();
	RUnivariateHelper(const Rcpp::Function& d, const Rcpp::Function& p,
		const Rcpp::Function& q, const Rcpp::Function& s);
	virtual ~RUnivariateHelper() { };

	virtual double d(double x, bool log = false) const;
	virtual double p(double q, bool lower = true, bool log = false) const;
	virtual double q(double p, bool lower = true, bool log = false) const;
	virtual bool s(double x) const;

private:
	Rcpp::Function _d;
	Rcpp::Function _p;
	Rcpp::Function _q;
	Rcpp::Function _s;
};

#endif
