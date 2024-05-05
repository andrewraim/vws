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
template <class T>
class FMMProposal
{
public:
	//' @description
	//' Constructor for FMMProposal.
	//' @param regions A list of objects whose class derives from \code{Region}.
	FMMProposal(std::vector<Region<T>> regions);

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
	typename std::set<Region<T>>::const_iterator get_regions() const;

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
	typename std::vector<T> r(unsigned int n = 1) const;

	typename std::pair<std::vector<T>, std::vector<unsigned int>> r_ext(unsigned int n = 1) const;

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

	void bifurcate(const Region<T>& r);

	void adapt(unsigned int N);

private:
	void recache();

	std::set<Region<T>> _regions; ////// TBD: Should these be sorted for subsequent searching?
	std::vector<double> _log_xi_upper;
	std::vector<double> _log_xi_lower;
	std::vector<bool> _bifurcatable;
};

// Recall that region volumes reflect where there mixture is further
// from the target: it takes into account both the weight difference
// and the probability of being in that region.
template <typename T>
void FMMProposal<T>::adapt(unsigned int N)
{
	for (unsigned int j = 0; j < N; j++) {
		unsigned int L = _regions.size();

		// Each region's contribution to the rejection rate
		Rcpp::NumericVector log_volume = rejection_bound_regions(true);

		std::vector<unsigned int> idx;
		for (unsigned int i = 0; i < L; i++) {
			if (_bifurcatable[i]) {
				idx.push_back(i);
			}
		}

		if (idx.size() == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		for (unsigned int l = 0; l < L; l++) {
			log_volume(l) *= std::pow(R_NegInf, _bifurcatable[l]);
		}

		unsigned int jdx = r_categ(1, log_volume, true);
		auto itr = _regions.begin();
		std::advance(itr, jdx);
		const Region<T>& r = *itr;

		// Split the target region and make another proposal with it
		const std::pair<std::unique_ptr<Region<T>>,std::unique_ptr<Region<T>>>& bif_out = r.bifurcate();
		typename std::set<Region<T>>::const_iterator itr_rm = _regions.find(r);
		_regions.erase(itr_rm);
		_regions.insert(*bif_out.left);
		_regions.insert(*bif_out.right);
		recache();
	}
}

template <class T>
FMMProposal<T>::FMMProposal(std::vector<Region<T>> regions)
	: _regions(), _log_xi_upper(), _log_xi_lower(), _bifurcatable()
{
	_regions.insert(regions.begin(), regions.end());
	recache();
}

template <class T>
void FMMProposal<T>::recache()
{
	_log_xi_upper.resize(_regions.size());
	_log_xi_lower.resize(_regions.size());
	_bifurcatable.resize(_regions.size());

	for (unsigned int j = 0; j < _regions.size(); j++) {
		_log_xi_upper[j] = _regions[j].xi_upper(true);
		_log_xi_lower[j] = _regions[j].xi_lower(true);
		_bifurcatable[j] = _regions[j].is_bifurcatable();
	}
}


template <class T>
std::vector<double>::const_iterator FMMProposal<T>::get_log_xi_upper() const
{
	return _log_xi_upper.begin();
}

template <class T>
std::vector<double>::const_iterator FMMProposal<T>::get_log_xi_lower() const
{
	return _log_xi_lower.begin();
}

template <class T>
std::vector<bool>::const_iterator FMMProposal<T>::get_bifurcatable() const
{
	return _bifurcatable.begin();
}

template <class T>
typename std::set<Region<T>>::const_iterator FMMProposal<T>::get_regions() const
{
	return _regions.begin();
}

template <class T>
Rcpp::NumericVector FMMProposal<T>::rejection_bound(bool log) const
{
	// Each region's contribution to the rejection rate bound
	Rcpp::NumericVector lxl(_log_xi_lower.begin(), _log_xi_lower.end());
	Rcpp::NumericVector lxu(_log_xi_upper.begin(), _log_xi_upper.end());
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);

	// Overall rejection rate bound
	out = vws::log_sum_exp(out);
	if (log) { return out; } else { return exp(out); }
}

template <class T>
Rcpp::NumericVector FMMProposal<T>::rejection_bound_regions(bool log) const
{
	// Each region's contribution to the rejection rate bound
	Rcpp::NumericVector lxl(_log_xi_lower.begin(), _log_xi_lower.end());
	Rcpp::NumericVector lxu(_log_xi_upper.begin(), _log_xi_upper.end());
	const Rcpp::NumericVector& out = log_sub2_exp(lxu, lxl) - log_sum_exp(lxu);
	if (log) { return out; } else { return exp(out); }
}

template <class T>
double FMMProposal<T>::nc(bool log) const
{
	Rcpp::NumericVector lxu(_log_xi_upper.begin(), _log_xi_upper.end());
	double out = vws::log_sum_exp(lxu);
	return log ? out : exp(out);
}

template <class T>
typename std::vector<T> FMMProposal<T>::r(unsigned int n) const
{
	return r_ext(n).first;
}

template <class T>
typename std::pair<std::vector<T>, std::vector<unsigned int>>
FMMProposal<T>::r_ext(unsigned int n) const
{
	unsigned int N = _regions.length();

	// Draw from the mixing weights, which are given on the log scale and
	// not normalized.
	Rcpp::NumericVector lp(_log_xi_upper.begin(), _log_xi_upper.end());
	const Rcpp::IntegerVector& idx = r_categ(n, lp, true);

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		unsigned int j = idx[i];
		x.push_back(_regions[j].r(1));
	}

	return std::make_pair(x, idx);
}

template <class T>
double FMMProposal<T>::d(const T& x, bool normalize, bool log) const
{
	double log_nc = 0;
	if (normalize) {
		log_nc = nc(true);
	}

	double out = NA_REAL;

	unsigned int N = _regions.size();

	// An idea
	// 1. Create a singleton region with x only.
	// 2. Find the region in r _regions with lower bound x.
	// 3. Make sure r contains x.
	//
	// We need a way to construct singleton regions of type Region<T> without
	// accessing the specific subtype.
	/*
	const Region<T>& x_rej _regions.first()->singleton(x);
	const Region<T>& reg = _regions.lower_bound(x_rej);
	if (!reg.s(x)) {
		Rcpp::stop("!reg.s(x)");
	}
	out = reg.w_major(x, true) + reg.d_base(x, true) - log_nc;
	*/

	// This search could be more efficient, but would need to be done in a
	// way that can support any kind of region. For example, if we can
	// define a "<" operator for region objects, we could consider a binary
	// search.
	/*
	for (unsigned int j = 0; j < N; j++) {
		const Region<T>& reg = _regions[j];
		if (reg.s(x)) {
			out = reg.w_major(x, true) + reg.d_base(x, true) - log_nc;
		}
	}
	*/

	return log ? out : exp(out);
}

template <class T>
double FMMProposal<T>::d_target_unnorm(const T& x, bool log) const
{
	const Region<T>& reg = (*_regions.begin());
	double out = reg.w(x, true) + reg.d_base(x, true);
	return log ? out : exp(out);
}

template <class T>
Rcpp::DataFrame FMMProposal<T>::summary() const
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

template <class T>
void FMMProposal<T>::print(unsigned int n) const
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
