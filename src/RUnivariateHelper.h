#ifndef R_UNIVARIATE_HELPER_H
#define R_UNIVARIATE_HELPER_H

#include "vws.h"

class RUnivariateHelper : public vws::UnivariateHelper<double>
{
public:
	RUnivariateHelper(const Rcpp::Function& pdf, const Rcpp::Function& cdf,
	 	const Rcpp::Function& quantile, const Rcpp::Function& supp);
	virtual ~RUnivariateHelper() { };

	virtual double pdf(double x, bool log = false) const;
	virtual double cdf(double x, bool lower = true, bool log = false) const;
	virtual double quantile(double x, bool lower = true, bool log = false) const;
	virtual bool supp(double x) const;

private:
	Rcpp::Function _pdf;
	Rcpp::Function _cdf;
	Rcpp::Function _quantile;
	Rcpp::Function _supp;
};

#endif
