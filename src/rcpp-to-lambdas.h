#ifndef RCPP_TO_LAMBDAS_H
#define RCPP_TO_LAMBDAS_H

#include "vws.h"

struct RejectionLambdas
{
	vws::uv_weight_function w;
	fntl::density d;
	fntl::cdf p;
	fntl::quantile q;
	vws::optimizer maxopt;
	vws::optimizer minopt;
};

RejectionLambdas rcpp_to_lambdas(const Rcpp::Function& w,
	const Rcpp::Function& d_base, const Rcpp::Function& p_base,
	const Rcpp::Function& q_base);

RejectionLambdas rcpp_to_lambdas(const Rcpp::Function& w,
	const Rcpp::Function& d_base, const Rcpp::Function& p_base,
	const Rcpp::Function& q_base, const Rcpp::Function& maxopt,
	const Rcpp::Function& minopt);

#endif
