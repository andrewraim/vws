#ifndef REJECTION_TUNE_H
#define REJECTION_TUNE_H

#include <RcppArmadillo.h>
#include "ConstSAEMajorizer.h"
#include "VWSStepOutput.h"
#include "local-util.h"

unsigned int tune(h, x, tol1, tol2)
{
	/*
	* Drop regions that contribute very little. But only if
	* overall bound is small enough, and dropping the region
	* does not put us back over tol1 threshold.
	*
	* Go through all mergeable pairs and see if we should merge
	* from our criteria. If we find any matches, do the merge and
	* start the process over. Keep doing this until we find no more
	* pairs to merge.
	*/

	bool repeat = true;
	unsigned int tunes = 0;

	while (repeat)
	{
		repeat = false;

		for (unsigned int j = 0; j < h.size(); j++)
		{
			const Rcpp::NumericVector& lbdd = h.bound_contrib(true);
			if (lbdd(j) >= log(tol2)) { continue; }

			const Rcpp::IntegerVector& idx = h.mergeable(j);

			for (unsigned int l = 0; l < idx.size(); l++)
			{
				// Merge regions j and l without modifying the proposal
				auto itr1 = std::next(h.regions_begin(), j);
				const R& merged = itr1->merge(*itr2);

				double lxu0 = xi_upper(true);
				double lxl0 = xi_lower(true);
				double lbdd0 = vws::log_sub2_exp(lxu0, lxl0) - vws::log_sum_exp(lxu0);

				// TBD: compute the overall bound if we replaced j and l with
				// the merged region.

				/*
				* If the bound of the merged region has not increased beyond
				* the threshold, carry out the same merge in the proposal.
				*/

				if (merged.log_bound < log(tol1)) {
					repeat = false;
					tunes++;
				}
			}
		}

		finished_merge = true
	}

	return tunes;
}

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
rejection_tune(FMMProposal<T,R>& h, unsigned int n,
	const rejection_args& args)
{

	rejection_result<T> out;

	unsigned int N_rejects = 0;
	bool accept = false;

	unsigned int max_rejects = args.max_rejects;
	unsigned int report = args.report;
	fntl::error_action action = args.action;
	double log_ratio_ub = std::exp(args.ratio_ub);

	// The constant M in the acceptance ratio is always M = 1.
	double log_M = 0;

	for (unsigned int i = 0; i < n && N_rejects <= max_rejects; i++)
	{
		accept = false;
		out.rejects.push_back(0L);

		while (!accept && N_rejects <= max_rejects)
		{
			double v = ::R::runif(0, 1);
			const T& x = h.r();
			double log_fx = h.d_target_unnorm(x);
			double log_hx = h.d(x, false, true);
			double log_ratio = log_fx - log_hx - log_M;

			if (log_ratio > log_ratio_ub) {
				Rcpp::stop("log_ratio %g exceeded %g; with x = %g, log f(x) = %g, and log h(x) = %g",
					log_ratio, log_ratio_ub, x, log_fx, log_hx);
			} else if (log(v) < log_ratio) {
				// Accept x as a draw from f(x)
				out.draws.push_back(x);
				accept = true;
			} else {
				// Reject x and refine the proposal
				N_rejects++;
				out.rejects[i]++;
			}

			// Report progress after `report` candidates
			unsigned int N_accepts = i + accept;
			if ((N_rejects + N_accepts) % report == 0) {
				Rprintf("%s - After %d candidates, %d accepts and %d rejects\n",
					timestamp().c_str(), N_accepts + N_rejects, N_accepts, N_rejects);
			}

			/* Consider tuning the proposal if the draw was rejected */
			if (!accept && h.bound(true) < log(tol1)) {
				// Identify regions that contribute very little and merge them
				// with larger regions.
				out.tunes(i) += tune_merge(h, x, tol1, tol2);
			} else if (!accept) {
				// Add the rejected draw as a knot
				h.({x});
				out.tunes(i)++;
			}
		}
	}

	if (N_rejects > max_rejects) {
		switch(action) {
			case fntl::error_action::STOP:
				Rcpp::stop("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			case fntl::error_action::WARNING:
				Rcpp::warning("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			case fntl::error_action::MESSAGE:
				Rprintf("Exceeded maximum number of rejects: %d\n", max_rejects);
				break;
			default:
				break;
		}
	}

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
inline rejection_result<T> rejection_tune(const FMMProposal<T,R>& h, unsigned int n = 1)
{
	rejection_args ctrl;
	return rejection_tune(h, n, ctrl);
}

#endif
