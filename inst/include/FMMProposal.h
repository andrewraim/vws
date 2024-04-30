#ifndef FMM_PROPOSAL_H
#define FMM_PROPOSAL_H

#include <Rcpp.h>
#include "log-sum-exp.h"
#include "categ.h"
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
	FMMProposal(ConstIterator itr);

	//' @description
	//' Access the vector \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	Rcpp::NumericVector get_xi_upper(bool log = false) const;

	//' @description
	//' Access the vector \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	Rcpp::NumericVector get_xi_lower(bool log = false) const;

	//' @description
	//' Access the vector \code{bifurcatable}.
	Rcpp::LogicalVector get_bifurcatable() const;

	//' @description
	//' Access the vector \code{regions}.
	ConstIterator get_regions() const;

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

	void bifurcate(const Region<T>& r);

	void adapt(unsigned int N);

private:
	void recache();

	std::set<Region<T>> _regions; ////// TBD: Should these be sorted for subsequent searching?
	Rcpp::NumericVector _log_xi_upper;
	Rcpp::NumericVector _log_xi_lower;
	Rcpp::LogicalVector _bifurcatable;
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

		const Rcpp::LogicalVector& is_bifurcatable = get_bifurcatable();

		const Rcpp::IntegerVector& idx = which(is_bifurcatable);
		if (idx.length() == 0) {
			Rcpp::warning("No regions left to bifurcate");
			break;
		}

		for (unsigned int l = 0; l < L; l++) {
			log_volume(l) *= std::pow(R_NegBin, is_bifurcatable);
		}

		unsigned int jdx = r_categ(1, log_volume, true);
		const Region& r = _regions[jdx];

		// Split the target region and make another proposal with it
		std::pair<Region<T>,Region<T>> bif_out = r.bifurcate();
		std::set<Region<T>>::iterator itr = _regions.find(r);
		_regions.erase(itr);
		_regions.insert(bif_out.left);
		_regions.insert(bif_out.right);
		recache();
	}
}

template <class T>
FMMProposal<T>::FMMProposal(ConstIterator itr)
	: _regions(), _log_xi_upper(), _log_xi_lower(), _bifurcatable()
{
	_regions.insert(itr.begin(), itr.end());
	recache();
}

template <class T>
void FMMProposal<T>::recache()
{
	_log_xi_upper.clear();
	_log_xi_lower.clear();
	_bifurcatable.clear();

	for (unsigned int j = 0; j < _regions.size(); j++) {
		_log_xi_upper[j] = _regions[j].xi_upper(true);
		_log_xi_lower[j] = _regions[j].xi_lower(true);
		_bifurcatable[j] = _regions[j].is_bifurcatable();
	}
}


template <class T>
Rcpp::NumericVector FMMProposal<T>::get_xi_upper(bool log) const
{
	const Rcpp::NumericVector& out = _log_xi_upper;
	if (log) { return out; } else { return exp(out); }
}

template <class T>
Rcpp::NumericVector FMMProposal<T>::get_xi_lower(bool log) const
{
	const Rcpp::NumericVector& out = _log_xi_lower;
	if (log) { return out; } else { return exp(out); }
}

template <class T>
Rcpp::LogicalVector FMMProposal<T>::get_bifurcatable() const
{
	return _bifurcatable;
}

template <class T>
std::vector<T> FMMProposal<T>::get_regions() const
{
	return _regions;
}

template <class T>
Rcpp::NumericVector FMMProposal<T>::rejection_bound(bool log) const
{
	// Each region's contribution to the rejection rate bound
	const Rcpp::NumericVector& out =
		log_sub2_exp(_log_xi_upper, _log_xi_lower) -
		log_sum_exp(_log_xi_upper);

	// Overall rejection rate bound
	out = vws::log_sum_exp(out);
	if (log) { return out; } else { return exp(out); }
}

template <class T>
Rcpp::NumericVector FMMProposal<T>::rejection_bound_regions(bool log) const
{
	// Each region's contribution to the rejection rate bound
	const Rcpp::NumericVector& out =
		log_sub2_exp(_log_xi_upper, _log_xi_lower) -
		log_sum_exp(_log_xi_upper);
	if (log) { return out; } else { return exp(out); }
}

template <class T>
double FMMProposal<T>::nc(bool log) const
{
	double out = vws::log_sum_exp(_log_xi_upper);
	return log ? out : exp(out);
}

template <class T>
std::vector<T> FMMProposal<T>::r(unsigned int n) const
{
	return r_ext(n).first;
}

template <class T>
std::pair<std::vector<T>, std::vector<unsigned int>>
FMMProposal<T>::r_ext(unsigned int n = 1) const
{
	unsigned int N = _regions.length();

	// Draw from the mixing weights, which are given on the log scale and
	// not normalized.
	idx = r_categ(n, _log_xi_upper, true);

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		unsigned int j = idx[i];
		x.push_back(_regions[j].r(1));
	}

	return std::pair<std::vector<T>, std::vector<unsigned int>>(x, idx);
}

template <class T>
double FMMProposal<T>::d(const T& x, bool normalize, bool log) const
{
	double log_nc = 0;
	if (normalize) {
		log_nc = nc(true);
	}

	double out = Rcpp::NA_REAL;

	unsigned int N = _regions.size();

	// This search could be more efficient, but would need to be done in a
	// way that can support any kind of region. For example, if we can
	// define a "<" operator for region objects, we could consider a binary
	// search.
	for (unsigned int j = 0; j < N; j++) {
		const Region<T>& reg = _regions[j];
		if (reg.s(x)) {
			out = reg.w_major(x, true) + reg.d_base(x, true) - log_nc;
		}
	}

	return log ? out : exp(out);
}

template <class T>
double FMMProposal<T>::d_target_unnorm(const T& x, bool log) const
{
	const Region<T>& reg = _regions[0];
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
	Rcpp::print(Rcpp::head(tbl, n));

	if (N > n) {
		printf("There are %d more regions not displayed\n", N - n);
	}
}

}

#endif
