// [[Rcpp::depends(vws, RcppFunctionalUtilities)]]
#include "vws.h"
#include "NormalHelper.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu,
	double sigma2, double lambda2, unsigned int N = 10,
	unsigned int max_rejects = 10000, unsigned int report_period = 1000)
{
	vws::RejectionControl ctrl;
	ctrl.max_rejects = max_rejects;
	ctrl.report_period = report_period;
	ctrl.max_rejects_action = vws::ErrorAction::STOP;

	const vws::weight_function& w =
    [&](double x, bool log = true) {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - mu, 2.0) / (2*sigma2);
		}
		return log ? out : std::exp(out);
	};

	NormalHelper helper(z, lambda2);
	vws::UnivariateConstRegion supp(0.0, R_PosInf, w, helper);
	vws::FMMProposal<double, vws::UnivariateConstRegion> h({ supp });

	h.adapt(N - 1);
	h.print(5);

	const vws::RejectionResult<double>& out = vws::rejection(h, n, ctrl);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects
	);
}
