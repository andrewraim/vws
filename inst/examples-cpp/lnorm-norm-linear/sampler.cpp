// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "normal-truncated.h"
#include "CustomLinearRegion.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report_period = report_period;
	args.max_rejects_action = vws::ErrorAction::STOP;

	double x = exp(mu - sigma2 + 1);
	CustomLinearRegion r1(1e-6, x, mu, sigma2, z, lambda2);
	CustomLinearRegion r2(x, 1e6, mu, sigma2, z, lambda2);
	const std::vector<CustomLinearRegion>& regions = { r1, r2 };

	vws::FMMProposal<double, CustomLinearRegion> h(regions);

	h.adapt(N - 1);
	h.print(100);

	const vws::rejection_result<double>& out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects
	);
}
