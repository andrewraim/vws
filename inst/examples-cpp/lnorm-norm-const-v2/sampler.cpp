// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "MyHelper.h"
#include "CustomConstRegion.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report_period = report_period;
	args.max_rejects_action = vws::error_action::STOP;

	MyHelper helper(z, lambda2);

	const vws::dfdb& w =
    [&](double x, bool log = true) {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - mu, 2.0) / (2*sigma2);
		}
		return log ? out : std::exp(out);
	};

    CustomConstRegion supp(0.0, R_PosInf, mu, sigma2, w, helper);

	const std::vector<CustomConstRegion>& regions = { supp };

	vws::FMMProposal<double, CustomConstRegion> h(regions);

	h.refine(N - 1);
	h.print(5);

	const vws::rejection_result<double>& out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects
	);
}
