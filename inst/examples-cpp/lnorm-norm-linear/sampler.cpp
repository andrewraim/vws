// [[Rcpp::depends(vws)]]
#include "vws.h"
//#include "CustomLinearRegion.h"
#include "normal_truncated.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	/*
	MyHelper helper(mu, sigma2, z, lambda2);
	CustomLinearRegion supp(0.0, R_PosInf, mu, sigma2, z, lambda2);

	const std::vector<CustomLinearRegion>& regions = { supp };

	vws::FMMProposal<double, CustomLinearRegion> h(regions);

	h.adapt(N - 1);
	h.print(5);

	vws::RejectionControl control(max_rejects, report_period, vws::MaxRejectsAction::stop);
	const std::pair<std::vector<double>, std::vector<unsigned int>>& out =
		vws::rejection(h, n, control);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejects") = out.second
	);
	*/

	return Rcpp::List::create();
}
