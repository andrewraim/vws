#ifndef VWS_FMM_PROPOSAL_H
#define VWS_FMM_PROPOSAL_H

#include <Rcpp.h>
#include "typedefs.h"
#include "log-sum-exp.h"
#include "categ.h"
#include "timestamp.h"
#include "Region.h"
#include <iterator>

namespace vws {

/*
* FMMProposal is a class which represents a VWS proposal: a finite mixture that
* implements certain operations. It takes two template arguments.
*
* - `T`: the data type of the support. For example, for a continuous
*   univariate support, we may use `double`.
* - `R`: a subclass of `Region` toused by the proposal. `Region` subclasses
*   implement any problem-specific logic used in the proposal.
*/
template <class T, class R>
class FMMProposal
{
public:
	/*
	* Constructor for FMMProposal based on a vector of one or more Regions.
	*/
	FMMProposal(const std::vector<R>& regions);

	/*
	* Constructor for FMMProposal based on a single Region. This is for
	* convenience in this commonly used special case.
	*/
	FMMProposal(const R& region);

	/*
	* TBD: copy constructor
	*/
	FMMProposal(const FMMProposal& p);

	/*
	* Access region characteristics
	*
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*
	* Return vectors which are ordered so that elements correspond to regions
	* $1, \ldots, N$.
	*
 	* - `xi_upper`: get $(\overline{\xi}_1, \ldots \overline{\xi}_N)$.
	* - `xi_lower`: get $(\underline{\xi}_1, \ldots \underline{\xi}_N)$.
	* - `bifurcatable`: vector of indicators for whether regions are
	*   bifurcatable.
	* - `pi`: mixing weights $(\pi_1, \ldots, \pi_N)$.
	* - `bound_contrib`: vector of contributions to the rejection
	*   bound. These sum to the value of `bound`.
	*/
	Rcpp::NumericVector xi_upper(bool log = true) const;
	Rcpp::NumericVector xi_lower(bool log = true) const;
	Rcpp::LogicalVector bifurcatable() const;
	Rcpp::NumericVector pi(bool log = false) const;
	Rcpp::NumericVector bound_contrib(bool log = false) const;

	// Return a vector of indices for regions which are mergeable to region i.
	// Does not include i itself and may be empty.
	Rcpp::IntegerVector mergeable(unsigned int i) const;

	/*
	* Access number of regions.
	*/
	unsigned int size() const {
		return _regions.size();
	}

	/*
	* Iterator to regions.
	*/
	typename std::set<R>::const_iterator begin() const {
		return _regions.begin();
	}
	typename std::set<R>::const_iterator end() const {
		return _regions.end();
	}

	// Rcpp::NumericVector::const_iterator log_xi_upper_begin() const {
	// 	return _log_xi_upper->begin();
	// }
	// Rcpp::NumericVector::const_iterator log_xi_upper_end() const {
	// 	return _log_xi_upper->end();
	// }
	// Rcpp::NumericVector::const_iterator log_xi_lower_begin() const {
	// 	return _log_xi_lower->begin();
	// }
	// Rcpp::NumericVector::const_iterator log_xi_lower_end() const {
	// 	return _log_xi_lower->end();
	// }
	// Rcpp::LogicalVector::const_iterator bifurcatable_begin() const {
	// 	return _bifurcatable->begin();
	// }
	// Rcpp::LogicalVector::const_iterator bifurcatable_end() const {
	// 	return _bifurcatable->end();
	// }

	/*
	* Upper bound for rejection probability.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double bound(bool log = false) const;

	/*
	* Normalizing constant for the proposal distribution.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double nc(bool log = false) const;

	/*
	*  Generate draws from proposal distribution.
	*  - `n`: desired number of draws.
	*  Returns a vector of `n` draws, or a single draw when `n` is unspecified.
	*/
	T r() const;
	std::vector<T> r(unsigned int n) const;

	/*
	*  Generate draws from proposal distribution.
	*  - `n`: desired number of draws.
	*  Returns a pair. First element is vector of draws. Second element is
	*  vector with indices, each represents the region from which the draw was
	*  taken.
	*/
	std::pair<std::vector<T>, std::vector<unsigned int>> r_ext(unsigned int n = 1) const;

	/*
	* Compute density of proposal distribution.
	* - `x`: argument of density.
	* - `normalize`: whether to compute normalized or unnormalized density.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double d(const T& x, bool normalize = true, bool log = false) const;

	/*
	* Compute $\overline{w}_(x)$.
	*
	* - `x`: a density argument.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double w_major(const T& x, bool log = true) const;

	/*
	* Compute $f_0(x) = w(x) g(x)$.
	*
	* - `x`: a density argument.
	* - `log`: if `true`, return value on the log-scale. Otherwise, return it
	*   on the original scale.
	*/
	double d_target_unnorm(const T& x, bool log = true) const;

	/*
	* Prepare summary table of regions which compose the proposal distribution.
	*/
	Rcpp::DataFrame summary() const;

	/*
	* Print summary table of regions which compose the proposal distribution.
	* Limit the display to `n` regions.
	*/
	void print(unsigned int n = 5) const;

	/*
	* TBD
	*/
	double merge(unsigned int i, unsigned int j, bool log = true);

	/*
	* Refine the proposal.
	*
	* - `knots`: vector of $N$ knots at which to bifurcate.
	*
	* This function modifies the proposal object in place.
	*
	* It returns a vector of $N+1$ values that represent rejection bounds
	* corresponding to the bifurcations from `knots`: the element with index
	* `j` is the bound after applying the first $j$ `knots`. The values are
	* returned on the log-scale if `log = true`.
	*/
	Rcpp::NumericVector refine(const std::vector<T>& knots, bool log = true);

	/*
	* Refine the proposal.
	*
	* - `N`: number of additional desired bifurcations.
	* - `tol`: tolerance for rejection bound.
	* - `greedy`: whether to use greedy selection to select regions for
	*   bifurcation.
	* - `report`: specifies the period in which progress should be reported
	*   (printed to the screen as a log message).
	*
	* This function modifies the proposal object in place.
	*
	* If `greedy` is `true`, the selected region to bifurcate will be the one
	* with the largest contribution to the rejection bound. Otherwise, the
	* selection is random, with probability proportional to regions'
	* contribution. Note that in either case, non-bifurcatable regions are
	* excluded from the candidates.
	*
	* If tolerance `tol` is achieved for the rejection bound before $N$
	* bifurcation steps are attempted, the operation is considered to be
	* complete.
	*
	* This function returns a vector of up to $N+1$ values that represent
	* rejection bounds (on the log-scale) that correspond to the bifurcation
	* steps: the element with index `j` corresponds to the bound at the $j$th
	* step. The values are returned on the log-scale if `log = true`.
	*/
	Rcpp::NumericVector refine(unsigned int N, double tol = 0,
		bool greedy = false, unsigned int report = fntl::uint_max, bool log = true);

protected:

	/*
	* Recompute internal auxiliary data structures.
	*
	* The following data structures are based on the set `_regions`:
	* - `_log_xi_upper`
	* - `_log_xi_lower`
	* - `_bifurcatable`
	*
	* Here we construct them from the current state of `_regions`.
	*/
	// void recache();

	std::set<R> _regions;
	// std::unique_ptr<Rcpp::NumericVector> _log_xi_upper;
	// std::unique_ptr<Rcpp::NumericVector> _log_xi_lower;
	// std::unique_ptr<Rcpp::LogicalVector> _bifurcatable;
};

template <class T, class R>
double FMMProposal<T,R>::merge(unsigned int i, unsigned int j, bool log)
{
	// Merge the regions with the given indices. Remove the two individual
	// regions and add the merged region. Call recache to update internal data
	// structures.
	auto itr1 = std::next(_regions.begin(), i);
	auto itr2 = std::next(_regions.begin(), j);
	const R& merged = itr1->merge(*itr2);
	_regions.erase(itr1);
	_regions.erase(itr2);
	_regions.insert(merged);
	// recache();

	double out = bound(true);
	return log ? out : std::exp(out);
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::refine(const std::vector<T>& knots, bool log)
{
	// unsigned int N = knots.size() + 1;

	std::vector<double> log_bdd_hist;

	log_bdd_hist.push_back(bound(true));

	for (unsigned int j = 0; j < knots.size(); j++)
	{
		// unsigned int L = _regions.size();

		const T& x = knots[j];
		const R& reg0 = _regions.begin()->singleton(x);

		typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(reg0);
		--itr_lower;

		if (itr_lower == _regions.end()) {
			Rcpp::warning("Could not find region which is upper bound for knots[%d]", j);
			continue;
		}

		// In this case, x was smaller (according to operator<) than all of the
		// regions in the set. So print a warning and continue.
		if (!itr_lower->s(x)) {
			Rcpp::warning("Could not find a region that contains knots[%d]", j);
			continue;
		}

		if (!itr_lower->is_bifurcatable()) {
			Rcpp::warning("Region that contains knots[%d] is not bifurcatable", j);
			continue;
		}

		const R& r = *itr_lower;

		// Bifurcate the selected region and replace the selected region with
		// the two new ones.
		const std::pair<R,R>& bif_out = r.bifurcate(x);
		typename std::set<R>::iterator itr = _regions.find(r);
		_regions.erase(itr);
		_regions.insert(bif_out.first);
		_regions.insert(bif_out.second);
		// recache();

		log_bdd_hist.push_back(bound(true));
	}

	Rcpp::NumericVector out(log_bdd_hist.begin(), log_bdd_hist.end());
	if (log) { return out; } else { return Rcpp::exp(out); }
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::refine(unsigned int N, double tol,
	bool greedy, unsigned int report, bool log)
{
	if (tol < 0) {
		Rcpp::stop("tol must be nonnegative");
	}

	std::vector<double> lbdd_hist;
	lbdd_hist.push_back(bound(true));

	for (unsigned int j = 0; j < N; j++) {
		// If we can beat the tolerance before we reach N steps, return now
		if (lbdd_hist[j] <= std::log(tol)) {
			break;
		}

		unsigned int L = _regions.size();

		// Each region's contribution to the rejection rate
		// const Rcpp::NumericVector& log_volume0 = bound_contrib(true);
		Rcpp::NumericVector log_volume(size());
		double lxusum = log_sum_exp(xi_upper(true));

		// Identify the regions which are bifurcatable; for the rest, we set
		// their log-volume to -inf so they will not be selected.
		unsigned int n_bif = 0;
		unsigned int l = 0;

		// for (unsigned int l = 0; l < L; l++) {
		for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
			double lxu = itr->xi_upper(true);
			double lxl = itr->xi_lower(true);
			double lrho = log_sub2_exp(lxu, lxl) - lxusum;

			if (std::isnan(lrho)) {
				Rcpp::stop("nan found in log_volume(%d)", l);
			}

			log_volume(l) = itr->is_bifurcatable() ? lrho : R_NegInf;
			n_bif += (log_volume(l) > R_NegInf);
			l++;
		}

		if (n_bif == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		unsigned int jdx;
		if (greedy) {
			jdx = Rcpp::which_max(log_volume);
		} else {
			jdx = r_categ(log_volume, true);
		}

		const R& r = *std::next(_regions.begin(), jdx);

		// Split the target region and make another proposal with it.
		//
		// TBD: we may not need to recache each time through the loop.
		// Perhaps we can insert into the end and zero out the removed entries?
		const std::pair<R,R>& bif_out = r.bifurcate();
		auto itr = _regions.find(r);
		_regions.erase(itr);
		_regions.insert(bif_out.first);
		_regions.insert(bif_out.second);
		// recache();

		lbdd_hist.push_back(bound(true));

		if (j % report == 0 && report < fntl::uint_max) {
			Rprintf("%s - After %d steps log Pr{rejection} <= %g\n",
				timestamp().c_str(), j, lbdd_hist[j+1]);
		}
	}

	Rcpp::NumericVector out(lbdd_hist.begin(), lbdd_hist.end());
	if (log) { return out; } else { return Rcpp::exp(out); }
}

template <class T, class R>
FMMProposal<T,R>::FMMProposal(const std::vector<R>& regions)
: _regions()
{
	_regions.insert(regions.begin(), regions.end());
	// recache();
}

template <class T, class R>
FMMProposal<T,R>::FMMProposal(const R& region)
: _regions()
{
	_regions.insert(region);
	// recache();
}

template <class T, class R>
FMMProposal<T,R>::FMMProposal(const FMMProposal& p)
: _regions()
{
	_regions.insert(p._regions.begin(), p._regions.end());
	// recache();
}

/*
template <class T, class R>
void FMMProposal<T,R>::recache()
{
	// TBD: This might be a good place to check for no overlaps, and perhaps to
	// make sure there are no gaps?

	unsigned int N = _regions.size();
	_log_xi_upper = std::make_unique<Rcpp::NumericVector>(N);
	_log_xi_lower = std::make_unique<Rcpp::NumericVector>(N);
	_bifurcatable = std::make_unique<Rcpp::LogicalVector>(N);

	unsigned int j = 0;
	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		(*_log_xi_upper)[j] = itr->xi_upper(true);
		(*_log_xi_lower)[j] = itr->xi_lower(true);
		(*_bifurcatable)[j] = itr->is_bifurcatable();
		j++;

		if (itr->xi_lower(true) > itr->xi_upper(true)) {
			Rcpp::stop("xi_lower > xi_upper");
		}
	}
}
*/

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::xi_upper(bool log) const
{
	unsigned int i = 0;
	Rcpp::NumericVector out(size());
	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		out(i) = itr->xi_upper(log);
		i++;
	}

	if (log) { return(out); } else { return Rcpp::exp(out); }
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::xi_lower(bool log) const
{
	unsigned int i = 0;
	Rcpp::NumericVector out(size());
	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		out(i) = itr->xi_lower(log);
		i++;
	}

	if (log) { return(out); } else { return Rcpp::exp(out); }
}

template <class T, class R>
Rcpp::LogicalVector FMMProposal<T,R>::bifurcatable() const
{
	unsigned int i = 0;
	Rcpp::LogicalVector out(size());
	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		out(i) = itr->is_bifurcateable();
		i++;
	}

	return out;
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::pi(bool log) const
{
	const Rcpp::NumericVector& lxu = xi_upper(true);
	const Rcpp::NumericVector& out = lxu - vws::log_sum_exp(lxu);
	if (log) { return(out); } else { return Rcpp::exp(out); }
}

template <class T, class R>
double FMMProposal<T,R>::bound(bool log) const
{
	// Overall rejection rate bound
	double out = log_sum_exp(bound_contrib(true));
	return log ? out : exp(out);
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::bound_contrib(bool log) const
{
	// Each region's contribution to the rejection rate bound.
	const Rcpp::NumericVector& lxl = xi_upper(true);
	const Rcpp::NumericVector& lxu = xi_lower(true);
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);
	if (log) { return(out); } else { return Rcpp::exp(out); }
}

template <class T, class R>
Rcpp::IntegerVector FMMProposal<T,R>::mergeable(unsigned int i) const
{
	unsigned int N = size();
	std::vector<unsigned int> out;

	auto itr1 = std::next(_regions.begin(), i);

	for (unsigned int j = 0; j < N; j++) {
		auto itr2 = std::next(_regions.begin(), j);
		if (i != j && itr1->is_mergeable(*itr2)) {
			out.push_back(j);
		}
	}

	return Rcpp::IntegerVector(out.begin(), out.end());
}

template <class T, class R>
double FMMProposal<T,R>::nc(bool log) const
{
	const Rcpp::NumericVector& lxu = xi_upper(true);
	double out = log_sum_exp(lxu);
	return log ? out : exp(out);
}

template <class T, class R>
T FMMProposal<T,R>::r() const
{
	return r_ext(1).first[0];
}

template <class T, class R>
typename std::vector<T> FMMProposal<T,R>::r(unsigned int n) const
{
	return r_ext(n).first;
}

template <class T, class R>
typename std::pair<std::vector<T>, std::vector<unsigned int>>
FMMProposal<T,R>::r_ext(unsigned int n) const
{
	// Draw from the mixing weights, which are given on the log scale and not
	// normalized.
	const Rcpp::NumericVector& lxu = xi_upper(true);
	const Rcpp::IntegerVector& idx = r_categ(n, lxu, true);

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		// const std::vector<T>& draws = _regions_vec[idx(i)].r(1);
		auto itr = std::next(_regions.begin(), idx(i));
		const std::vector<T>& draws = itr->r(1);
		x.push_back(draws[0]);
	}

	return std::make_pair(x, Rcpp::as<std::vector<unsigned int>>(idx));
}

template <class T, class R>
double FMMProposal<T,R>::d(const T& x, bool normalize, bool log) const
{
	double lnc = normalize ? nc(true) : 0;

	/*
	* Search for the region containing x using std::set. This is a bit more
	* complicated than just iterating over all the regions until we find the
	* right one, but it should take something like O(log N) time rather than
	* O(N), where N is the number of regions.
	*
	* 1. Create a singleton region `x_singleton` from the point x.
	* 2. Find the region r in  _regions whose upper bound is x. This should be
	*    one position past the region we want.
	* 3. Rewind the iterator one step. This should yield the region with x.
	* 4. Make sure this region contains x; otherwise, throw an error.
	*/
	const R& x_singleton = _regions.begin()->singleton(x);
	typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(x_singleton);
	--itr_lower;

	double out;

	if (itr_lower == _regions.end()) {
		// Could not find region which is upper bound for x. Assume that x is
		// outside of the support
		out = R_NegInf;
		// Rprintf("FMMProposal<T,R>::d, case 1: out = %g\n", out);
	} else if (!itr_lower->s(x)) {
		// Could not find a region that contains x. Assume that x is outside of
		// the support.
		out = R_NegInf;
		// Rprintf("FMMProposal<T,R>::d, case 2: out = %g\n", out);
	} else {
		out = itr_lower->w_major(x, true) + itr_lower->d_base(x, true) - lnc;
		// Rprintf("FMMProposal<T,R>::d, case 3: out = %g  log-wmajor = %g  log-base-density = %g  log-nc = %g\n",
		// 	out, itr_lower->w_major(x, true), itr_lower->d_base(x, true), lnc);
	}

	return log ? out : exp(out);
}

template <class T, class R>
double FMMProposal<T,R>::d_target_unnorm(const T& x, bool log) const
{
	auto itr = _regions.begin();
	double out = itr->w(x, true) + itr->d_base(x, true);
	return log ? out : exp(out);
}

template <class T, class R>
double FMMProposal<T,R>::w_major(const T& x, bool log) const
{
	/*
	* Search for the region containing x using std::set, similar to the `d`
	* method.
	*/
	const R& x_singleton = _regions.begin()->singleton(x);

	typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(x_singleton);
	--itr_lower;

	double out;

	if (itr_lower == _regions.end()) {
		// Could not find region which is upper bound for x. Assume that x is
		// outside of the support
		out = R_NegInf;
	} else if (!itr_lower->s(x)) {
		// Could not find a region that contains x. Assume that x is outside of
		// the support.
		out = R_NegInf;
	} else {
		out = itr_lower->w_major(x, true);
	}

	return log ? out : exp(out);
}

template <class T, class R>
Rcpp::DataFrame FMMProposal<T,R>::summary() const
{
	unsigned int N = _regions.size();
	Rcpp::StringVector v1(N);
	Rcpp::NumericVector v2(N);
	Rcpp::NumericVector v3(N);
	Rcpp::NumericVector v4(N);
	double lden = log_sum_exp(xi_upper(true));

	unsigned int j = 0;

	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		double lxl = itr->xi_lower(true);
		double lxu = itr->xi_upper(true);
		v1[j] = itr->description();
		v2[j] = lxu;
		v3[j] = lxl;
		v4[j] = log_sub2_exp(lxu, lxl) - lden;
		j++;
	}

	return Rcpp::DataFrame::create(
		Rcpp::Named("Region") = v1,
		Rcpp::Named("log_xi_upper") = v2,
		Rcpp::Named("log_xi_lower") = v3,
		Rcpp::Named("log_volume") = v4
	);
}

template <class T, class R>
void FMMProposal<T,R>::print(unsigned int n) const
{
	const Rcpp::DataFrame& tbl = summary();
	unsigned int N = tbl.nrows();

	const Rcpp::IntegerVector& idx = Rcpp::seq(0, std::min(n,N) - 1);

	const Rcpp::StringVector& v1 = tbl["Region"];
	const Rcpp::NumericVector& v2 = tbl["log_xi_upper"];
	const Rcpp::NumericVector& v3 = tbl["log_xi_lower"];
	const Rcpp::NumericVector& v4 = tbl["log_volume"];

	const Rcpp::DataFrame& tbl_head = Rcpp::DataFrame::create(
		Rcpp::Named("Region") = v1[idx],
		Rcpp::Named("log_xi_upper") = v2[idx],
		Rcpp::Named("log_xi_lower") = v3[idx],
		Rcpp::Named("log_volume") = v4[idx]
	);

	Rprintf("FMM Proposal with %d regions\n", N);

	Rcpp::print(tbl_head);

	if (N > n) {
		Rprintf("There are %d more regions not displayed\n", N - n);
	}
}

}

#endif
