#include <Rcpp.h>
#include "RUnivariateHelper.h"
#include "UnivariateConstRegion.h"

//' @export
// [[Rcpp::export]]
double my_test_function(const Rcpp::Environment& e)
{
	// TBD: do we need to do anything else to check that the cast was appropriate?
	SEXP op = e[".pointer"];
	const RUnivariateHelper* p = static_cast<RUnivariateHelper*>(R_ExternalPtrAddr(op));
	return p->pdf(3, true);
}

RCPP_MODULE(VWSModule) {

Rcpp::class_<RUnivariateHelper>("RUnivariateHelper")
	.constructor<Rcpp::Function, Rcpp::Function,
		Rcpp::Function, Rcpp::Function>()
	.method("pdf", &RUnivariateHelper::pdf)
	.method("cdf", &RUnivariateHelper::cdf)
	.method("quantile", &RUnivariateHelper::quantile)
	.method("supp", &RUnivariateHelper::supp);

/*
Rcpp::class_<vws::UnivariateConstRegion>("UnivariateConstRegion")
	.constructor<double, double, Rcpp::Function, Rcpp::Environment>()
	.method("r", &vws::UnivariateConstRegion::r);
*/

}
