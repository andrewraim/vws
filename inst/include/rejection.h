#ifndef REJECTION_H
#define REJECTION_H

#include <Rcpp.h>
#include "logger.h"
#include "Region.h"
#include "RejectionControl.h"

namespace vws {

//' Vertical Weighted Strips Rejection Sampler
//'
//' Accept-reject algorithm using our proposal for weighted distributions.
//'
//' @param h An \code{fmm_proposal} object
//' @param n Number of desired draws
//' @param control A control object from \code{rejection_control}
//'
//' @return
//' A list whose structure depends on
//' \code{extra_outputs}. If \code{extra_outputs = FALSE}, the list is of length
//' \code{n} where each element represents one draw. If
//' \code{extra_outputs = TRUE}, the list contains the following named elements:
//'
//' \item{draws}{is list is of length \code{n} where each element represents one
//' draw.}
//' \item{rejects}{is a vector of counts. The \eqn{i}th element is the number of
//' rejections before the \eqn{i}th successful draw.}
//' \item{h}{returns the \code{h} that was input to the sampler.}
//'
//' Note that \code{h} was originally intended to show adaptation during
//' sampling, but that is currently not done.
//'
//' @examples
//' # Define base distribution and weight function
//' g = normal_univariate_helper(0, 5)
//' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
//'
//' # Set up support
//' support = UnivariateConstRegion$new(-Inf, 10, w, g)
//' regions = support$bifurcate()
//'
//' # Create a finite mixture proposal
//' h = FMMProposal$new(regions)
//' h$rejection_bound()
//' h$rejection_bound(byregion = TRUE)
//'
//' out = rejection(h, n = 1000)
//' print(out |> unlist())
//'
//' out = rejection(h, n = 1000, rejection_control(extra_outputs = TRUE))
//' print(out$draws |> unlist())
//' print(out$rejects)
//'
//' @name rejection
//' @export
template <typename T, typename R>
std::pair<std::vector<T>, std::vector<unsigned int>>
rejection(const FMMProposal<T,R>& h, unsigned int n, const RejectionControl& control)
{
	std::vector<T> out;
	std::vector<unsigned int> rejects(n, 0L);

	unsigned int N_rejects = 0;
	bool accept = false;

	unsigned int max_rejects = control.get_max_rejects();
	unsigned int report_period = control.get_report_period();
	ErrorAction max_rejects_action = control.get_max_rejects_action();

	// The constant M in the acceptance ratio is always M = 1.
	double log_M = 0;

	for (unsigned int i = 0; i < n; i++) {
		accept = false;
		while (!accept && N_rejects < max_rejects) {
			double v = ::R::runif(0, 1);
			const std::vector<T>& draws = h.r(1);
			const T& x = draws[0];
			double log_fx = h.d_target_unnorm(x);
			double log_hx = h.d(x, false, true);
			double log_ratio = log_fx - log_hx - log_M;

			if (log(v) < log_ratio) {
				// Accept x as a draw from f(x)
				out.push_back(x);
				accept = true;
			} else {
				// Reject x and adapt the proposal
				N_rejects++;
				rejects[i]++;
			}

			// Report progress after `report` candidates
			unsigned int N_accepts = i + accept;
			if ((N_rejects + N_accepts) % report_period == 0) {
				logger("After %d candidates, %d accepts and %d rejects\n",
					N_accepts + N_rejects, N_accepts, N_rejects);
			}
		}
	}

	if (N_rejects >= max_rejects) {
		switch(max_rejects_action) {
			case ErrorAction::STOP:
		    	Rcpp::stop("Reached maximum number of rejects: %d\n", max_rejects);
		    	break;
		    case ErrorAction::WARNING:
		    	Rcpp::warning("Reached maximum number of rejects: %d\n", max_rejects);
		    	break;
		    case ErrorAction::MESSAGE:
		    	Rprintf("Reached maximum number of rejects: %d\n", max_rejects);
		    	break;
		}
	}

	return std::make_pair(out, rejects);
}

template <typename T, typename R>
std::pair<std::vector<T>, std::vector<unsigned int>>
rejection(const FMMProposal<T,R>& h, unsigned int n)
{
	// Set defaults for control object
	// - Up to 1000 rejects per requested draw.
	// - Do not report progress.
	// - Halt if we reach max_rejects before accepting the requested n draws.
	unsigned int max_rejects = 1000 * n;
	unsigned int report_period = n + 1;
	ErrorAction max_rejects_action = ErrorAction::STOP;

	vws::RejectionControl control(max_rejects, report_period, max_rejects_action);
	return rejection(h, n, control);
}

}

#endif
