#include <Rcpp.h>
#include "RUnivariateHelper.h"

RCPP_MODULE(VWSModule) {

Rcpp::class_<RUnivariateHelper>("RUnivariateHelper")
	.constructor()
	.constructor<Rcpp::Function, Rcpp::Function, Rcpp::Function, Rcpp::Function>()
	.method("d", &RUnivariateHelper::d)
	.method("p", &RUnivariateHelper::p)
	.method("q", &RUnivariateHelper::q)
	.method("s", &RUnivariateHelper::s);

}
