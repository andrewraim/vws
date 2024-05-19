// [[Rcpp::depends(vws, RcppFunctionalUtilities)]]
#include "vws.h"
#include "NormalHelper.h"

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	vws::RejectionControl control(max_rejects, report_period, vws::ErrorAction::STOP);

	const vws::weight_function& w =
    [&](double x, bool log = true) {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - mu, 2.0) / (2*sigma2);
		}
		return log ? out : exp(out);
	};

	NormalHelper helper(z, lambda2);
	vws::UnivariateConstRegion supp(0.0, R_PosInf, w, helper);
	vws::FMMProposal<double, vws::UnivariateConstRegion> h({ supp });

	h.adapt(N - 1);
	h.print(5);

	// const std::pair<std::vector<double>, std::vector<unsigned int>>& out =
	const auto& out = vws::rejection(h, n, control);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejects") = out.second
	);
}
