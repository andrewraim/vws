#ifndef ADAPT_H
#define ADAPT_H

#include <Rcpp.h>
#include "log-sum-exp.h"
#include "which.h"
#include "FMMProposal.h"

namespace vws {

//' Adapt Using Midpoint
//'
//' Adapt an FMM proposal using a midpoint rule-of-thumb.
//'
//' @param h An object of class \code{FMMProposal}.
//' @param N Number of additional mixture components after adaptation.
//' @param report_period Report progress each time this many candidates are proposed.
//'
//' @examples
//' g = normal_univariate_helper(0, 5)
//' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
//' support = UnivariateConstRegion$new(-Inf, 10, w, g)
//' regions = list(support)
//'
//' # Adapt the proposal
//' h_init = FMMProposal$new(regions)
//' adapt_out = adapt(h_init, N = 100)
//' h = adapt_out$h
//'
//' # Create a finite mixture proposal
//' h$rejection_bound()
//' h$rejection_bound(byregion = TRUE)
//'
//' out = rejection(h, n = 200)
//' print(out |> unlist())
//'
//' out = rejection(h, n = 200, rejection_control(extra_outputs = TRUE))
//' print(out$draws |> unlist())
//' print(out$rejects)
//'
//' @name adapt
//' @export
template <typename T>
std::tuple<FMMProposal<T>, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>
adapt(const FMMProposal<T>& h_init, unsigned int N, unsigned int report_period)
{
	FMMProposal<T> h = h_init;

	Rcpp::NumericVector log_ub_hist(N+1);
	Rcpp::NumericVector log_lb_hist(N+1);
	Rcpp::NumericVector log_bdd_hist(N+1);

	log_ub_hist(0) = log_sum_exp(h.get_xi_upper(true));
	log_lb_hist(0) = log_sum_exp(h.get_xi_lower(true));
	log_bdd_hist(0) = h.rejection_bound(true);

	if (report_period <= N) {
		logger("Initial log Pr{rejection} <= %g\n", log_bdd_hist(0));
	}

	for (unsigned int j = 0; j < N; j++) {
		// Recall that region volumes reflect where there mixture is further
		// from the target: it takes into account both the weight difference
		// and the probability of being in that region.
		unsigned int L = h.get_regions().length();
		const Rcpp::NumericVector& log_xi_upper = h.get_xi_upper(true);
		const Rcpp::NumericVector& log_xi_lower = h.get_xi_lower(true);

		// Each region's contribution to the rejection rate
		Rcpp::NumericVector log_volume = h.rejection_bound_regions(true);

		const Rcpp::LogicalVector& is_bifurcatable = h.get_bifurcatable();

		const Rcpp::IntegerVector& idx = which(is_bifurcatable);
		if (idx.length() == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		for (unsigned int i = 0; i < L; i++) {
			log_volume(i) *= std::pow(R_NegBin, is_bifurcatable);
		}

		unsigned int jdx = r_categ(1, log_volume, true);
		const Region& reg = h.get_regions()[jdx];

		// Split the target region and make another proposal with it
		h.bifurcate(reg);

		log_ub_hist(j+1) = log_sum_exp(h.get_xi_upper(true));
		log_lb_hist(j+1) = log_sum_exp(h.get_xi_lower(true));
		log_bdd_hist(j+1) = h.rejection_bound(true);

		if (j % report_period == 0) {
			logger("After %d steps log Pr{rejection} <= %g\n", j, log_bdd_hist(j+1));
		}
	}

	return std::make_tuple(h, log_ub_hist, log_lb_hist, log_bdd_hist);
}

}

#endif
