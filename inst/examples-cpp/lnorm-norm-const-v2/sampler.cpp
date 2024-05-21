// [[Rcpp::depends(vws, RcppFunctionalUtilities)]]
#include "vws.h"
#include "MyHelper.h"
#include "CustomConstRegion.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	vws::RejectionControl ctrl;
	ctrl.max_rejects = max_rejects;
	ctrl.report_period = report_period;
	ctrl.max_rejects_action = vws::ErrorAction::STOP;

	MyHelper helper(z, lambda2);

	const vws::weight_function& w =
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

	h.adapt(N - 1);
	h.print(5);

	const RejectionResult& out = vws::rejection(h, n, ctrl);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects
	);
}
