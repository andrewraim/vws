#ifndef VWS_FMM_PROPOSAL_H
#define VWS_FMM_PROPOSAL_H

#include <Rcpp.h>
#include "typedefs.h"
#include "log-sum-exp.h"
#include "categ.h"
#include "which.h"
#include "logger.h"
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
	Rcpp::NumericVector::const_iterator get_log_xi_upper() const;

	//' @description
	//' Access the vector \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	Rcpp::NumericVector::const_iterator get_log_xi_lower() const;

	//' @description
	//' Access the vector \code{bifurcatable}.
	Rcpp::LogicalVector::const_iterator get_bifurcatable() const;

	//' @description
	//' Access the vector \code{regions}.
	typename std::set<R>::const_iterator get_regions() const;

	//' @description
	//' Upper bound for rejection probability.
	//' @param byregion If \code{TRUE}, compute bound by region. Otherwise compute
	//' total.
	//' @param log If \code{TRUE} compute result on log-scale.
	//' @return A vector of size \code{N} or a single scalar.
	double rejection_bound(bool log = false) const;

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

	void adapt(unsigned int N, double tol = R_NegInf, unsigned int report = uint_max);

private:
	void recache();

	std::set<R> _regions;
	std::vector<R> _regions_vec;
	std::unique_ptr<Rcpp::NumericVector> _log_xi_upper;
	std::unique_ptr<Rcpp::NumericVector> _log_xi_lower;
	std::unique_ptr<Rcpp::LogicalVector> _bifurcatable;
};

// Recall that region volumes reflect where there mixture is further from the
// target: it takes into account both the weight difference and the probability
// of being in that region.
template <class T, class R>
void FMMProposal<T,R>::adapt(unsigned int N, double tol, unsigned int report)
{
	Rcpp::NumericVector log_bdd_hist(N+1);
	log_bdd_hist(0) = rejection_bound(true);

	for (unsigned int j = 0; j < N; j++) {
		// If we can beat the tolerance before we reach N steps, return now
		if (log_bdd_hist(j) <= tol) {
			return;
		}

		unsigned int L = _regions.size();

		// Each region's contribution to the rejection rate
		Rcpp::NumericVector log_volume = rejection_bound_regions(true);

		// printf("Checkpoint: log_volume before transform\n");
		// Rcpp::print(log_volume);

		unsigned int n_bif = 0;
		for (unsigned int l = 0; l < L; l++) {
			log_volume(l) = (*_bifurcatable)(l) ? log_volume(l) : R_NegInf;
			n_bif += (*_bifurcatable)(l);
		}

		if (n_bif == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		// printf("Checkpoint: log_volume after transform\n");
		// Rcpp::print(log_volume);

		unsigned int jdx = r_categ(log_volume, true);
		// printf("Checkpoint: jdx = %d\n", jdx);
		const R& r = _regions_vec[jdx];

		// Split the target region and make another proposal with it.
		//
		// TBD: we may not need to recache each time through the loop.
		// Perhaps we can insert into the end and zero out the removed entries?
		// printf("Checkpoint: begin split\n");
		const std::pair<R,R>& bif_out = r.bifurcate();
		// printf("Checkpoint: bifurcate\n");
		typename std::set<R>::iterator itr = _regions.find(r);
		// printf("Checkpoint: erase\n");
		_regions.erase(itr);
		// printf("Checkpoint: insert first\n");
		_regions.insert(bif_out.first);
		// printf("Checkpoint: insert second\n");
		_regions.insert(bif_out.second);
		// printf("Checkpoint: recache\n");
		recache();

		log_bdd_hist(j+1) = rejection_bound(true);

		// printf("Checkpoint: end split\n");
		if (j % report == 0 && report < uint_max) {
			logger("After %d steps log Pr{rejection} <= %g\n", j, log_bdd_hist(j+1));
		}
	}
}

template <class T, class R>
FMMProposal<T,R>::FMMProposal(const std::vector<R>& regions)
: _regions(), _regions_vec()
{
	_regions.insert(regions.begin(), regions.end());
	recache();
}

template <class T, class R>
void FMMProposal<T,R>::recache()
{
	// TBD: This might be a good place to check for no overlaps, and perhaps to
	// make sure there are no gaps?

	unsigned int N = _regions.size();
	_regions_vec.clear();
	_log_xi_upper = std::make_unique<Rcpp::NumericVector>(_regions.size());
	_log_xi_lower = std::make_unique<Rcpp::NumericVector>(_regions.size());
	_bifurcatable = std::make_unique<Rcpp::LogicalVector>(_regions.size());

	unsigned int j = 0;
	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		_regions_vec.push_back(*itr);
		(*_log_xi_upper)[j] = itr->get_xi_upper(true);
		(*_log_xi_lower)[j] = itr->get_xi_lower(true);
		(*_bifurcatable)[j] = itr->is_bifurcatable();
		j++;
	}
}

template <class T, class R>
Rcpp::NumericVector::const_iterator FMMProposal<T,R>::get_log_xi_upper() const
{
	return _log_xi_upper->begin();
}

template <class T, class R>
Rcpp::NumericVector::const_iterator FMMProposal<T,R>::get_log_xi_lower() const
{
	return _log_xi_lower->begin();
}

template <class T, class R>
Rcpp::LogicalVector::const_iterator FMMProposal<T,R>::get_bifurcatable() const
{
	return _bifurcatable->begin();
}

template <class T, class R>
typename std::set<R>::const_iterator FMMProposal<T,R>::get_regions() const
{
	return _regions.begin();
}

template <class T, class R>
double FMMProposal<T,R>::rejection_bound(bool log) const
{
	// Each region's contribution to the rejection rate bound
	const Rcpp::NumericVector& lxl = *_log_xi_lower;
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	const Rcpp::NumericVector& log_bound = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);

	// Overall rejection rate bound
	double out = log_sum_exp(log_bound);
	if (log) { return out; } else { return exp(out); }
}

template <class T, class R>
Rcpp::NumericVector FMMProposal<T,R>::rejection_bound_regions(bool log) const
{
	// Each region's contribution to the rejection rate bound.
	const Rcpp::NumericVector& lxl = *_log_xi_lower;
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);
	if (log) { return out; } else { return exp(out); }
}

template <class T, class R>
double FMMProposal<T,R>::nc(bool log) const
{
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	double out = log_sum_exp(lxu);
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
	unsigned int N = _regions.size();

	// Draw from the mixing weights, which are given on the log scale and
	// not normalized.
	const Rcpp::NumericVector& lxu = *_log_xi_upper;
	const Rcpp::IntegerVector& idx = r_categ(n, lxu, true);

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		const std::vector<T>& draws = _regions_vec[idx(i)].r(1);
		x.push_back(draws[0]);
	}

	return std::make_pair(x, Rcpp::as<std::vector<unsigned int>>(idx));
}

template <class T, class R>
double FMMProposal<T,R>::d(const T& x, bool normalize, bool log) const
{
	double log_nc = normalize ? nc(true) : 0;

	/*
	* Search for the region containing x using std::set. This should require
	* something like O(log N) time, where N is the number of regions.
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

	if (itr_lower == _regions.end()) {
		Rcpp::stop("Could not find region with point x");
	}

	if (!itr_lower->s(x)) {
		x_singleton.print();
		itr_lower->print();
		Rcpp::stop("!itr_lower->s(x)");
	}

	double out = itr_lower->w_major(x, true) + itr_lower->d_base(x, true) - log_nc;

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
	unsigned int N = _regions.size();
	Rcpp::StringVector v1(N);
	Rcpp::NumericVector v2(N);
	Rcpp::NumericVector v3(N);

	unsigned int j = 0;

	for (auto itr = _regions.begin(); itr != _regions.end(); ++itr) {
		v1[j] = itr->description();
		v2[j] = itr->get_xi_upper(true);
		v3[j] = itr->get_xi_lower(true);
		j++;
	}

	return Rcpp::DataFrame::create(
		Rcpp::Named("Region") = v1,
		Rcpp::Named("log_xi_upper") = v2,
		Rcpp::Named("log_xi_lower") = v3
	);
}

template <class T, class R>
void FMMProposal<T,R>::print(unsigned int n) const
{
	// unsigned int N = _regions.size();
	const Rcpp::DataFrame& tbl = summary();
	unsigned int N = tbl.nrows();

	const Rcpp::IntegerVector& idx = Rcpp::seq(0, std::min(n,N) - 1);

	const Rcpp::StringVector& v1 = tbl["Region"];
	const Rcpp::NumericVector& v2 = tbl["log_xi_upper"];
	const Rcpp::NumericVector& v3 = tbl["log_xi_lower"];

	const Rcpp::DataFrame& tbl_head = Rcpp::DataFrame::create(
		Rcpp::Named("Region") = v1[idx],
		Rcpp::Named("log_xi_upper") = v2[idx],
		Rcpp::Named("log_xi_lower") = v3[idx]
	);

	printf("FMM Proposal with %d regions\n", N);

	Rcpp::print(tbl_head);

	if (N > n) {
		printf("There are %d more regions not displayed\n", N - n);
	}
}

}

#endif
