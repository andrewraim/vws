#ifndef ADAPT_H
#define ADAPT_H

#include "FMMProposal.h"
#include "log-sum-exp.h"

namespace vws {

//' Adapt Using Midpoint
//'
//' Adapt an FMM proposal using a midpoint rule-of-thumb.
//'
//' @param h An object of class \code{FMMProposal}.
//' @param N Number of additional mixture components after adaptation.
//' @param report Report progress each time this many candidates are proposed.
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
std::tuple<FMMProposal, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>
adapt(const FMMProposal<T>& h_init, unsigned int N, unsigned int report = N+1)
{
	FMMProposal<T> h = h_init;

	log_ub_hist = Rcpp::NumericVector(N+1);
	log_lb_hist = Rcpp::NumericVector(N+1);
	log_bdd_hist = Rcpp::NumericVector(N+1);

	log_ub_hist(0) = log_sum_exp(h.get_xi_upper(true));
	log_lb_hist(0) = log_sum_exp(h.get_xi_lower(true));
	log_bdd_hist(0) = h.rejection_bound(true);

	if (report <= N) {
		logger("Initial log Pr{rejection} <= %g\n", log_bdd_hist(0));
	}

	for (unsigned int j = 0; j < N; j++) {
		// Recall that region volumes reflect where there mixture is further
		// from the target: it takes into account both the weight difference
		// and the probability of being in that region.
		L = h.get_regions().length();
		const Rcpp::NumericVector& log_xi_upper = h.xi_upper(true);
		const Rcpp::NumericVector& log_xi_lower = h.xi_lower(true);

		// Each region's contribution to the rejection rate
		double log_volume = h.rejection_bound(true, true);

		idx = Rcpp::which(h.get_bifurcatable());
		if (idx.length() == 0) {
			warning("No regions left to bifurcate");
			break;
		}

		jdx = r_categ(log_volume(idx), true);
		reg = h.get_regions()[[idx[jdx]]];

		// Split the target region and make another proposal with it
		// TBD: implement this split
		bif_out = reg.bifurcate();
		regions_new = append(h.get_regions()[-idx[jdx]], bif_out);
		h = FMMProposal(regions_new);

		log_ub_hist(j+1) = log_sum_exp(h.xi_upper(true));
		log_lb_hist(j+1) = log_sum_exp(h.xi_lower(true));
		log_bdd_hist(j+1) = h.rejection_bound(true);

		if (j %% report == 0) {
			logger("After %d steps log Pr{rejection} <= %g\n", j, log_bdd_hist(j+1));
		}
	}

	return std::make_tuple(h, log_ub_hist, log_lb_hist, log_bdd_hist);
}

}

#endif
