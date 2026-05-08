#ifndef REJECTION_TUNE_H
#define REJECTION_TUNE_H

namespace vws {

/*
* For regions that contribute very little (i.e., whose log-bound contribution
* is smaller than tol2), merge them with other regions. Only do this merge if
* the overall bound is small enough, and dropping the region would not put us
* back over tol1 threshold.
*
* Return `true` if we complete a merge; the caller can then call again to
* attempt another merge if they like. This is to prevent mistakes in
* iterators and indices when the underlying data is modified.
*/
template <typename T, typename R>
unsigned int tune_merge_one(FMMProposal<T,R>& h, const T& x, double tol1, double tol2)
{
	/*
	* Go through all mergeable pairs and check if we should merge (using our
	* criteria) without modifying the proposal. If we find a match, do the
	* merge in the proposal and return.
	*/

	for (unsigned int j = 0; j < h.size(); j++)
	{
		const Rcpp::NumericVector& lbdd = h.bound_contrib(true);
		if (lbdd(j) >= log(tol2)) { continue; }

		const Rcpp::IntegerVector& idx = h.mergeable(j);
		auto itr1 = std::next(h.regions_begin(), j);

		for (unsigned int l = 0; l < idx.size(); l++)
		{
			// Merge regions j and l without modifying the proposal
			auto itr2 = std::next(h.regions_begin(), idx(j));
			const R& merged = itr1->merge(*itr2);

			// Compute the overall bound if we replaced regions j and l
			// with a merge of those two regions.
			double lbdd0 = vws::log_add2_exp(h.bound(true), merged.bound(true));
			lbdd0 = vws::log_sub2_exp(lbdd0, lbdd(j));
			lbdd0 = vws::log_sub2_exp(lbdd0, lbdd(l));

			// If the total bound with the merged region has not decreased
			// below threshold tol1, perform the merge in the proposal.
			if (lbdd0 >= log(tol1))
			{
				h.merge(j, l);
				return true;
			}
		}
	}

	return false;
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
rejection_tune(FMMProposal<T,R>& h, unsigned int n, const rejection_args& args)
{
	rejection_result<T> out;

	unsigned int N_rejects = 0;
	unsigned int N_tunes = 0;
	bool accept = false;

	unsigned int max_rejects = args.max_rejects;
	unsigned int report = args.report;
	fntl::error_action action = args.action;
	double log_ratio_ub = std::exp(args.ratio_ub);
	double tol1 = args.tol1;
	double tol2 = args.tol2;

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
				Rcpp::stop("log_ratio %g exceeded %g; with x = %g, "
					"log f(x) = %g, and log h(x) = %g",
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
				Rprintf("%s - %d candidates  %d accepts  %d rejects  %d tunes\n",
					timestamp().c_str(), N_accepts + N_rejects, N_accepts,
					N_rejects, N_tunes);
			}

			/*
			 *  Consider tuning the proposal.
			 *
			 * If the draw was rejected and the bound for rejection
			 * rate is "sufficient", see if we can merge some of the regions
			 * to simplify the proposal. We merge two regions at a time,
			 * repeating until the next merge would make the rejection bound
			 * "insufficient".
			 *
			 *
			 * Otherwise, if the draw was rejected and the bound for rejection
			 * rate is is "insufficient", partition using the rejected x.
			 *
			 * Note that we do not do both on a single rejection.
			*/

			if (!accept && h.bound(true) >= log(tol1)) {
				// Identify regions that contribute very little and merge them
				// with other regions.
				bool repeat = true;

				/*
				* Go through all mergeable pairs and check if we should merge, using our
				* criteria, without modifying the proposal. If we find a match, do the
				* merge in the proposal and start the process over. Repeat until we find no
				* more pairs to merge.
				*/
				while (repeat) {
					bool merged = tune_merge_one(h, x, tol1, tol2);
					repeat = merged;
					out.tunes[i] += merged;
				}
			} else if (!accept && h.bound(true) < log(tol1)) {
				// Add the rejected draw as a knot
				const std::vector<T>& knots = { x };
				h.refine(knots);
				out.tunes[i]++;
				N_tunes++;
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

}

#endif
