// [[Rcpp::depends(vws, RcppFunctionalUtilities)]]
#include "vws.h"
#include "MyHelper.h"
#include "CustomConstRegion.h"

double callFunction(const std::function<double(const Rcpp::NumericVector&)>& f)
{
	const Rcpp::NumericVector& x = Rcpp::NumericVector::create(1, 2, 3);
    return f(x);
}

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	printf("Checkpoint: Creating helper and supp\n");
	MyHelper helper(mu, sigma2, z, lambda2);
	CustomConstRegion supp(0.0, R_PosInf, mu, sigma2, helper);

	printf("Checkpoint: Creating regions\n");
	const std::vector<CustomConstRegion>& regions = { supp };

	printf("Checkpoint: Creating proposal\n");
	vws::FMMProposal<double, CustomConstRegion> h(regions);

	printf("Checkpoint: Adapting proposal\n");
	h.adapt(N - 1);
	h.print(5);

	printf("Checkpoint: running sampler\n");
	vws::RejectionControl control(max_rejects, report_period, vws::ErrorAction::STOP);
	const std::pair<std::vector<double>, std::vector<unsigned int>>& out =
		vws::rejection(h, n, control);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejects") = out.second
	);
}
