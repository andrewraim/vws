#ifndef FMM_PROPOSAL_H
#define FMM_PROPOSAL_H

#include <Rcpp.h>
#include "log-sum-exp.h"
#include "categ.h"
#include "which.h"
#include "Region.h"
#include <iterator>

namespace vws {

//' FMMProposal
//'
//' An R6 class which represents a VWS proposal: a finite mixture that some
//' specific certain operations.
//'
//' @details
//' \itemize{
//' \item The list \code{regions} represents
//' (\eqn{\mathcal{D}_1, \ldots, \mathcal{D}_N}).
//' Each entry should be an object based on an R6  subclass of \code{Region}.
//' \item The numeric vectors \code{log_xi_upper} and \code{log_xi_lower}
//' represent
//' \eqn{(\log \overline{\xi}_1, \ldots, \log \overline{\xi}_N)}
//' and
//' \eqn{(\log \underline{\xi}_1, \ldots, \log \underline{\xi}_N)},
//' respectively.
//' \item The logical vector \code{bifurcatable} indicates whether each
//' region can be bifurcated or not.
//' }
//'
//' @export
template <class T, class R>
class FMMProposal
{
public:
	//' @description
	//' Constructor for FMMProposal.
	//' @param regions A list of objects whose class derives from \code{Region}.
	FMMProposal(const std::vector<R>& regions);

	//' @description
	//' Access the vector \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	std::vector<double>::const_iterator get_log_xi_upper() const;

	//' @description
	//' Access the vector \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	std::vector<double>::const_iterator get_log_xi_lower() const;

	//' @description
	//' Access the vector \code{bifurcatable}.
	std::vector<bool>::const_iterator get_bifurcatable() const;

	//' @description
	//' Access the vector \code{regions}.
	typename std::set<R>::const_iterator get_regions() const;

	//' @description
	//' Upper bound for rejection probability.
	//' @param byregion If \code{TRUE}, compute bound by region. Otherwise compute
	//' total.
	//' @param log If \code{TRUE} compute result on log-scale.
	//' @return A vector of size \code{N} or a single scalar.
	Rcpp::NumericVector rejection_bound(bool log = false) const;

	Rcpp::NumericVector rejection_bound_regions(bool log = false) const;

	//' @description
	//' Normalizing constant for proposal distribution.
	//' @param log If \code{TRUE} compute result on log-scale
	double nc(bool log = false) const;

	//' @description
	//' Generate draws from proposal distribution.
	//' @param n Number of draws to generate.
	//' @param indices If \code{TRUE}, return indices of mixture components selected
	//' during draws.
	//' @return A list which each element is a saved draw.
	std::vector<T> r(unsigned int n = 1) const;

	std::pair<std::vector<T>, std::vector<unsigned int>> r_ext(unsigned int n = 1) const;

	//' @description
	//' Compute density of proposal distribution.
	//' @param x A list where each element is a density value to evaluate.
	//' @param normalize If \code{TRUE}, apply the normalizing constant; otherwise
	//' do not.
	//' @param log If \code{TRUE}, return density results on the log-scale.
	//' @return A vector or scalar of density values corresponding to \code{x}.
	double d(const T& x, bool normalize = true, bool log = false) const;

	//' @description
	//' Compute \eqn{f_0(x) = w(x) g(x)} for each element of the given \code{x}
	//' which is a list whose elements are density values.
	//' @param x a list where each element represents one argument.
	//' @param log If \code{TRUE}, return density results on the log-scale.
	//' @return A vector or scalar of density values corresponding to \code{x}.
	double d_target_unnorm(const T& x, bool log = true) const;

	//' @description
	//' Summary table of regions which compose the proposal distribution. Returns a
	//' data frame.
	Rcpp::DataFrame summary() const;

	//' @description
	//' Display summary table of regions which compose the proposal distribution.
	//' Limit the display to \code{n} regions.
	//' @param n Number of regions to print.
	void print(unsigned int n = 5) const;

	void bifurcate(const R& r);

	void adapt(unsigned int N);

private:
	void recache();

	std::set<R> _regions;
	std::vector<R> _regions_vec;
	std::vector<double> _log_xi_upper;
	std::vector<double> _log_xi_lower;
	std::vector<bool> _bifurcatable;
};

// Recall that region volumes reflect where there mixture is further
// from the target: it takes into account both the weight difference
// and the probability of being in that region.
template <class T, class R>
void FMMProposal<T,R>::adapt(unsigned int N)
{
	for (unsigned int j = 0; j < N; j++) {
		unsigned int L = _regions.size();

		// Each region's contribution to the rejection rate
		Rcpp::NumericVector log_volume = rejection_bound_regions(true);

		unsigned int n_bif = 0;
		for (unsigned int l = 0; l < L; l++) {
			// Rprintf("pre:  log_volume(%d) = %g\n", l, log_volume(l));
			// log_volume(l) *= std::pow(R_NegInf, _bifurcatable[l]);
			log_volume(l) = _bifurcatable[l] ? log_volume(l) : R_NegInf;
			n_bif += _bifurcatable[l];
			// Rprintf("post: log_volume(%d) = %g\n", l, log_volume(l));
		}

		if (n_bif == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		unsigned int jdx = r_categ(log_volume, true);
		// Rprintf("jdx = %d", jdx);
		const R& r = _regions_vec[jdx];

		// Split the target region and make another proposal with it
		// TBD: we may not need to recache... perhaps we can insert into the
		// end and zero out the removed entries? Especially for this function
		// in particular, we may not need to recache each time through the loop.
		const R& bif1 = r.bifurcate_first();
		const R& bif2 = r.bifurcate_second();
		typename std::set<R>::iterator itr = _regions.find(r);
		_regions.erase(itr);
		_regions.insert(bif1);
		_regions.insert(bif2);
		recache();
	}
}

template <class T, class R>
FMMProposal<T,R>::FMMProposal(const std::vector<R>& regions)
	: _regions(), _regions_vec(), _log_xi_upper(), _log_xi_lower(), _bifurcatable()
{
	_regions.insert(regions.begin(), regions.end());
	recache();
}

template <class T, class R>
void FMMProposal<T,R>::recache()
{
	// TBD: This might be a good place to check for no overlaps, and perhaps to
	// make sure there are no gaps?

	_regions_vec.clear();
	_log_xi_upper.clear();
	_log_xi_lower.clear();
	_bifurcatable.clear();

	typename std::set<R>::const_iterator itr = _regions.begin();
	for (; itr != _regions.end(); ++itr) {
		_regions_vec.push_back(*itr);
		_log_xi_upper.push_back(itr->get_xi_upper(true));
		_log_xi_lower.push_back(itr->get_xi_lower(true));
		_bifurcatable.push_back(itr->is_bifurcatable());
	}
}

template <class T, class R>
std::vector<double>::const_iterator FMMProposal<T,R>::get_log_xi_upper() const
{
	return _log_xi_upper.begin();
}

template <class T, class R>
std::vector<double>::const_iterator FMMProposal<T,R>::get_log_xi_lower() const
{
	return _log_xi_lower.begin();
}

template <class T, class R>
std::vector<bool>::const_iterator FMMProposal<T,R>::get_bifurcatable() const
{
	return _bifurcatable.begin();
}

template <class T, class R>
typename std::set<R>::const_iterator FMMProposal<T,R>::get_regions() const
{
	return _regions.begin();
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::rejection_bound(bool log) const
{
	// Each region's contribution to the rejection rate bound
	Rcpp::NumericVector lxl(_log_xi_lower.begin(), _log_xi_lower.end());
	Rcpp::NumericVector lxu(_log_xi_upper.begin(), _log_xi_upper.end());
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);

	// Overall rejection rate bound
	out = vws::log_sum_exp(out);
	if (log) { return out; } else { return exp(out); }
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::rejection_bound_regions(bool log) const
{
	// Each region's contribution to the rejection rate bound
	Rcpp::NumericVector lxl(_log_xi_lower.begin(), _log_xi_lower.end());
	Rcpp::NumericVector lxu(_log_xi_upper.begin(), _log_xi_upper.end());
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);

	// Rcpp::print(lxl);
	// Rcpp::print(lxu);
	// Rcpp::print(out);

	if (log) { return out; } else { return exp(out); }
}

template <class T, class R>
double FMMProposal<T,R>::nc(bool log) const
{
	Rcpp::NumericVector lxu(_log_xi_upper.begin(), _log_xi_upper.end());
	double out = vws::log_sum_exp(lxu);
	return log ? out : exp(out);
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
	// Rprintf("r_ext: Checkpoint 1\n");

	unsigned int N = _regions.size();

	// Rprintf("r_ext: Checkpoint 2\n");

	// Draw from the mixing weights, which are given on the log scale and
	// not normalized.
	Rcpp::NumericVector lp(_log_xi_upper.begin(), _log_xi_upper.end());
	const Rcpp::IntegerVector& idx = r_categ(n, lp, true);

	// Rcpp::print(lp);
	// Rcpp::print(idx);
	// Rcpp::stop("PAUSE!");

	// Rprintf("r_ext: Checkpoint 3\n");

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		// Rprintf("r_ext: Checkpoint 3.1\n");
		unsigned int j = idx(i);
		// Rprintf("r_ext: Checkpoint 3.2, j = %d\n", j);
		const std::vector<T>& draws = _regions_vec[j].r(1);
		// Rprintf("r_ext: Checkpoint 3.3\n");
		x.push_back(draws[0]);
		// Rprintf("r_ext: Checkpoint 3.4\n");
	}

	// Rprintf("r_ext: Checkpoint 4\n");

	return std::make_pair(x, Rcpp::as<std::vector<unsigned int>>(idx));
}

template <class T, class R>
double FMMProposal<T,R>::d(const T& x, bool normalize, bool log) const
{
	double log_nc = normalize ? nc(true) : 0;

	// Rprintf("d: Checkpoint 1\n");

	double out = NA_REAL;

	// Rprintf("d: Checkpoint 2\n");

	unsigned int N = _regions.size();

	// Rprintf("d: Checkpoint 3\n");

	// An idea
	// 1. Create a singleton region with x only.
	// 2. Find the region in r _regions with lower bound x.
	// 3. Make sure r contains x.
	//
	// Construct a singleton region with x to find which partition contains x.
	// To do this, use the upper_bound function in the set class, which returns
	// the region *after* the one we want. Then we rewind it one position to
	// get the target region.
	const R& x_rej = _regions.begin()->singleton(x);
	// Rprintf("d: Checkpoint 4, N = %d\n", N);
	// x_rej.print();
	// _regions.begin()->print();
	// Rprintf("*_regions.begin() < x_rej = %d\n", (*_regions.begin()) < x_rej);
	typename std::set<R>::const_iterator itr_lower = _regions.upper_bound(x_rej);
	--itr_lower;
	// Rprintf("d: Checkpoint 5\n");
	if (itr_lower == _regions.end()) {
		// itr_lower->print();
		Rcpp::stop("Could not find region with point x");
	}
	// Rprintf("d: Checkpoint 5.1\n");
	if (!itr_lower->s(x)) {
		Rcpp::stop("!itr_lower->s(x)");
	}
	// Rprintf("d: Checkpoint 6\n");
	out = itr_lower->w_major(x, true) + itr_lower->d_base(x, true) - log_nc;

	// Rprintf("d: Checkpoint 7\n");

	// This search could be more efficient, but would need to be done in a
	// way that can support any kind of region. For example, if we can
	// define a "<" operator for region objects, we could consider a binary
	// search.
	/*
	for (unsigned int j = 0; j < N; j++) {
		const R& reg = _regions[j];
		if (reg.s(x)) {
			out = reg.w_major(x, true) + reg.d_base(x, true) - log_nc;
		}
	}
	*/

	return log ? out : exp(out);
}

template <class T, class R>
double FMMProposal<T,R>::d_target_unnorm(const T& x, bool log) const
{
	const R& r = *_regions.begin();
	double out = r.w(x, true) + r.d_base(x, true);
	return log ? out : exp(out);
}

template <class T, class R>
Rcpp::DataFrame FMMProposal<T,R>::summary() const
{
	/*
	tbl = data.frame(
		Region = Map(function(x) { x$description() }, _regions) |> unlist(),
		log_xi_upper = Map(function(x) { x$xi_upper(log = TRUE) }, _regions) |> unlist(),
		log_xi_lower = Map(function(x) { x$xi_lower(log = TRUE) }, _regions) |> unlist()
	)
	return tbl;
	*/
	return Rcpp::DataFrame();
}

template <class T, class R>
void FMMProposal<T,R>::print(unsigned int n) const
{
	unsigned int N = _regions.size();
	printf("FMM Proposal with %d regions (display is unsorted)\n", N);

	const Rcpp::DataFrame& tbl = summary();
	// Rcpp::print(Rcpp::head(tbl, n));
	Rcpp::print(tbl);

	if (N > n) {
		printf("There are %d more regions not displayed\n", N - n);
	}
}

}

#endif
