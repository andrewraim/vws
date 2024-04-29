#ifndef FMMPROPOSAL_H
#define FMMPROPOSAL_H

#include <Rcpp.h>
#include "log-sum-exp.h"
#include "categ.h"

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
template <class T> class FMMProposal {
private:
	vector<Region> _regions; ////// TBD: Should these be sorted for subsequent searching?
	Rcpp::NumericVector _log_xi_upper;
	Rcpp::NumericVector _log_xi_lower;
	vector<bool> _bifurcatable;

public:
	//' @description
	//' Constructor for FMMProposal.
	//' @param regions A list of objects whose class derives from \code{Region}.
	FMMProposal(const std::vector<Region>& regions);

	//' @description
	//' Access the vector \eqn{\overline{\xi}_1, \ldots, \overline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	Rcpp::NumericVector xi_upper(bool log = false) const;

	//' @description
	//' Access the vector \eqn{\underline{\xi}_1, \ldots, \underline{\xi}_N}.
	//' @param log If \code{TRUE} compute result on log-scale.
	Rcpp::NumericVector xi_lower(bool log = false) const;

	//' @description
	//' Access the vector \code{bifurcatable}.
	bool get_bifurcatable() const;

	//' @description
	//' Access the vector \code{regions}.
	vector<T> get_regions() const;

	//' @description
	//' Upper bound for rejection probability.
	//' @param byregion If \code{TRUE}, compute bound by region. Otherwise compute
	//' total.
	//' @param log If \code{TRUE} compute result on log-scale.
	//' @return A vector of size \code{N} or a single scalar.
	Rcpp::NumericVector rejection_bound(bool byregion = false, bool log = false) const;

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

	std::pair<std::vector<T>, std::vector<unsigned int>>
	r_ext(unsigned int n = 1) const;

	//' @description
	//' Compute density of proposal distribution.
	//' @param x A list where each element is a density value to evaluate.
	//' @param normalize If \code{TRUE}, apply the normalizing constant; otherwise
	//' do not.
	//' @param log If \code{TRUE}, return density results on the log-scale.
	//' @return A vector or scalar of density values corresponding to \code{x}.
	double d(const T& x, normalize = true, log = false) const;

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
};


FMMProposal(const std::vector<Region>& regions)
: _regions(), _log_xi_upper(), _log_xi_lower(), _bifurcatable()
{
	_regions = regions;
	for (unsigned int j = 0; j < _regions.size(); j++) {
		log_xi_upper.push_back(_regions.xi_upper(true));
		log_xi_lower.push_back(_regions.xi_lower(true));
		_bifurcatable.push_back(_regions.is_bifurcatable());
	}
}

Rcpp::NumericVector FMMProposal::xi_upper(bool log = false) const
{
	out = _log_xi_upper;
	if (log) { return(out); } else { return(exp(out)); }
}


Rcpp::NumericVector FMMProposal::xi_lower(bool log = false) const
{
	out = _log_xi_lower;
	if (log) { return(out); } else { return(exp(out)); }
}

bool FMMProposal::get_bifurcatable() const
{
	return _bifurcatable;
}

vector<T> FMMProposal::get_regions() const
{
	return _regions;
}

Rcpp::NumericVector FMMProposal::rejection_bound(bool byregion = false, bool log = false) const
{
	// Each region's contribution to the rejection rate bound
	const Rcpp::NumericVector& out =
		log_sub2_exp(_log_xi_upper, _log_xi_lower) -
		log_sum_exp(_log_xi_upper);

	if (!byregion) {
		// Overall rejection rate bound
		out = vws::log_sum_exp(out);
	}

	return log ? out : exp(out);
}

double FMMProposal::nc(bool log = false) const
{
	out = vws::log_sum_exp(_log_xi_upper);
	return log ? out : exp(out);
}

std::vector<T> FMMProposal::r(unsigned int n = 1) const
{
	return r_ext(n).first;
}

std::pair<std::vector<T>, std::vector<unsigned int>>
FMMProposal::r_ext(unsigned int n = 1) const
{
	unsigned int N = length(_regions);

	// Draw from the mixing weights, which are given on the log scale and
	// not normalized.
	idx = r_categ(n, p = _log_xi_upper, log_p = true);

	// Draw the values from the respective mixture components.
	std::vector<T> x;
	for (unsigned int i = 0; i < n; i++) {
		j = idx[i];
		x.push_back(_regions[j].r(1));
	}

	return std::pair<std::vector<T>, std::vector<unsigned int>>(x, idx);
}

double FMMProposal::d(const T& x, normalize = true, log = false) const
{
	double log_nc = 0;
	if (normalize) {
		log_nc = nc(log = TRUE)
	}

	double out = Rcpp::NA_REAL;

	unsigned int N = _regions.size();

	// This search could be more efficient, but would need to be done in a
	// way that can support any kind of region. For example, if we can
	// define a "<" operator for region objects, we could consider a binary
	// search.
	for (unsigned int j = 0; j < N; j++) {
		reg = _regions[j];
		if (reg.s(x[i])) {
			out = reg.w_major(x, true) + reg.d_base(x, true) - log_nc;
		}
	}

	return log ? out : exp(out);
}

double FMMProposal::d_target_unnorm(const T& x, bool log = true) const
{
	reg = _regions[0];
	out = Map(function(z) { reg.w(z, true) + reg.d_base(z, true) }, x);
	if (log) { return(unlist(out)); } else { return(exp(unlist(out))); }
}

Rcpp::DataFrame FMMProposal::summary() const
{
	tbl = data.frame(
		Region = Map(function(x) { x$description() }, _regions) |> unlist(),
		log_xi_upper = Map(function(x) { x$xi_upper(log = TRUE) }, _regions) |> unlist(),
		log_xi_lower = Map(function(x) { x$xi_lower(log = TRUE) }, _regions) |> unlist()
	)
	return tbl;
}

void FMMProposal::print(unsigned int n = 5) const
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
