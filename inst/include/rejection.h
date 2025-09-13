#ifndef VWS_REJECTION_H
#define VWS_REJECTION_H

#include <Rcpp.h>
#include "logger.h"
#include "Region.h"
#include "FMMProposal.h"
#include "result.h"
#include "typedefs.h"
#include "rejection-args.h"

namespace vws {

/*
*  Accept-reject algorithm using a VWS proposal.
*
*  - `h`: a VWS proposal.
*  - `n`: number of desired draws.
*  - `args`: additional arguments
*
*  Returns a structure with saved draws and rejection counts.
*/
template <typename T, typename R>
inline rejection_result<T>
rejection(const FMMProposal<T,R>& h, unsigned int n, const rejection_args& args)
{
	std::vector<T> draws;
	std::vector<unsigned int> rejects;

	unsigned int N_rejects = 0;
	bool accept = false;

	// Rprintf("Rejection Checkpoint 1\n");

	unsigned int max_rejects = args.max_rejects;
	unsigned int report = args.report;
	error_action action = args.action;
	double log_ratio_ub = std::exp(args.ratio_ub);

	// Rprintf("Rejection checkpoint 2\n");

	// The constant M in the acceptance ratio is always M = 1.
	double log_M = 0;

	for (unsigned int i = 0; i < n && N_rejects <= max_rejects; i++)
	{
		accept = false;
		rejects.push_back(0L);

		while (!accept && N_rejects <= max_rejects)
		{
			// Rprintf("Rejection checkpoint 2.1\n");
			double v = ::R::runif(0, 1);
			// Rprintf("Rejection checkpoint 2.2\n");
			const T& x = h.r();
			// Rprintf("Rejection checkpoint 2.3\n");
			double log_fx = h.d_target_unnorm(x);
			// Rprintf("Rejection checkpoint 2.4\n");
			double log_hx = h.d(x, false, true);
			// Rprintf("Rejection checkpoint 2.5\n");
			double log_ratio = log_fx - log_hx - log_M;

			// Rprintf("Rejection checkpoint 2.2\n");
			if (log_ratio > log_ratio_ub) {
				Rcpp::stop("log_ratio %g exceeded %g; with x = %g, log f(x) = %g, and log h(x) = %g",
					log_ratio, log_ratio_ub, x, log_fx, log_hx);
			} else if (log(v) < log_ratio) {
				// Accept x as a draw from f(x)
				draws.push_back(x);
				accept = true;
			} else {
				// Reject x and refine the proposal
				N_rejects++;
				rejects[i]++;
			}

			// Rprintf("Rejection checkpoint 2.3\n");
			// Report progress after `report` candidates
			unsigned int N_accepts = i + accept;
			if ((N_rejects + N_accepts) % report == 0) {
				logger("After %d candidates, %d accepts and %d rejects\n",
					N_accepts + N_rejects, N_accepts, N_rejects);
			}
		}
	}

	// Rprintf("Rejection checkpoint 3\n");

	if (N_rejects > max_rejects) {
		switch(action) {
			case error_action::STOP:
				Rcpp::stop("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			case error_action::WARNING:
				Rcpp::warning("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			case error_action::MESSAGE:
				Rprintf("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			default:
				break;
		}
	}

	// Rprintf("Rejection checkpoint 4\n");

	rejection_result<T> out;
	out.draws = draws;
	out.rejects = rejects;

	// Rprintf("Rejection checkpoint 5\n");
	return out;
}

/*
*  Accept-reject algorithm using VWS proposal.
*
*  - `h`: a VWS proposal.
*  - `n`: number of desired draws.
*
*  Returns a structure with saved draws and rejection counts.
*/
template <typename T, typename R>
inline rejection_result<T> rejection(const FMMProposal<T,R>& h, unsigned int n)
{
	rejection_args ctrl;
	return rejection(h, n, ctrl);
}

}

#endif
