// [[Rcpp::depends(vws)]]
#include "vws.h"
#include "normal-truncated.h"
#include "CustomLinearRegion.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	// printf("Create regions\n");

	double x = exp(mu - sigma2 + 1);

	// printf("x = %g\n", x);

	CustomLinearRegion r1(1e-3, x, mu, sigma2, z, lambda2);
	CustomLinearRegion r2(x, 1e6, mu, sigma2, z, lambda2);

	const std::vector<CustomLinearRegion>& regions = { r1, r2 };

	// printf("Create proposal\n");

	vws::FMMProposal<double, CustomLinearRegion> h(regions);

	// printf("Adapt proposal\n");

	h.adapt(N - 1);
	h.print(100);

	// printf("Take a sample\n");

	// const std::vector<double>& draws = h.r(100);
	// Rcpp::NumericVector draws2(draws.begin(), draws.end());
	// Rcpp::print(draws2);

	// printf("Create proposal\n");

	vws::RejectionControl control(max_rejects, report_period, vws::MaxRejectsAction::stop);
	const std::pair<std::vector<double>, std::vector<unsigned int>>& out =
		vws::rejection(h, n, control);

	// printf("Packing up results\n");

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejects") = out.second
	);

	return Rcpp::List::create();
}
