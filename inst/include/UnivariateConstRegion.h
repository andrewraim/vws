#ifndef UNIVARIATE_CONST_REGION_H
#define UNIVARIATE_CONST_REGION_H

#include <Rcpp.h>
#include <memory>
#include <RcppFunctionalUtilities.h>
#include "Region.h"
#include "UnivariateHelper.h"

namespace vws {

//' Univariate Region with Constant Majorizer
//'
//' An R6 class which represents a region based on univariate intervals with a
//' constant majorizer for the weight function. This version is for continuous
//' supports.
//'
//' @field w Weight function for the target distribution. Its expected interface
//' is \code{w(x, log = TRUE)} so that results are returned on the log-scale by
//' default.
//'
//' @examples
//' # Define base distribution and weight function
//' g = normal_univariate_helper(mean = 0, sd = 5)
//' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
//'
//' reg = UnivariateConstRegion$new(-Inf, 10, w, g)
//' print(reg)
//'
//' out = reg$bifurcate(0)
//' print(out[[1]])
//' print(out[[2]])
//'
//' @export
class UnivariateConstRegion : public Region<double>
{
protected:
	double _a;
	double _b;
	const UnivariateHelper<double>* _helper;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;

public:
	//' @param a Lower limit of interval.
	//' @param b Upper limit of interval.
	//' @param w Weight function for the target distribution.
	//' @param g An object created by \code{univariate_helper}.
	UnivariateConstRegion(double a, double b, const UnivariateHelper<double>& helper);

	//' @description
	//' Density function \eqn{g} for the base distribution.
	//' @param x Density argument.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double d_base(const double& x, bool log = false) const;

	//' @description
	//' Generate a draw from \eqn{g_j} specific to this region.
	//' @param n Number of draws to generate.
	//' @return A list of draws, with one draw per list element.
	std::vector<double> r(unsigned int n) const;

	//' @description
	//' Density of \eqn{g_j} specific to this region.
	//' @param x Density argument.
	double d(const double& x, bool log = false) const;

	//' @description
	//' Test if given \code{x} is in the support for the \eqn{g_j} specific to
	//' this region.
	//' @param x Density argument.
	bool s(const double& x) const;

	double w(const double& x, bool log = true) const;

	//' @description
	//' Return a logical value indicating whether this region is bifurcatable.
	bool is_bifurcatable() const;

	//' @description
	//' The quantity \eqn{\overline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double get_xi_upper(bool log = true) const;

	//' @description
	//' The quantity \eqn{\underline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double get_xi_lower(bool log = true) const;

	//' @description
	//' A string that describes the region.
	std::string description() const;

	//' @description
	//' Print a description of the region.
	void print() const;

	//' @description
	//' Maximize or minimize the function \eqn{w(x)} over this region. Optimization
	//' is carried out with \code{optim} using arguments
	//' \code{method = "Nelder-Mead"} and
	//' \code{control = list(maxit = 100000, warn.1d.NelderMead = FALSE)}.
	//' @param maximize logical; if \code{TRUE} do maximization. Otherwise do
	//' minimization.
	//' @param log logical; if \code{TRUE} return optimized value of \eqn{\log w(x)}.
	//' Otherwise return optimized value of \eqn{w(x)}.
	double optimize(bool maximize = true, bool log = true) const;

	//' @description
	//' Majorized weight function \eqn{\overline{w}_j} for this region.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double w_major(const double& x, bool log = true) const;

	double midpoint() const;

	std::pair<UnivariateConstRegion,UnivariateConstRegion> bifurcate() const;

	//' @description
	//' Bifurcate this region into two regions. Use \code{x} as the bifurcation
	//' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
	//' @param x An optional bifurcation point.
	std::pair<UnivariateConstRegion,UnivariateConstRegion> bifurcate(const double& x) const;

	UnivariateConstRegion singleton(const double& x) const;

	bool operator<(const UnivariateConstRegion& x) const {
		return _a < x._a;
	}

	bool operator==(const UnivariateConstRegion& x) const {
		return _a == x._a && _b == x._b;
	}

	const UnivariateConstRegion& operator=(const UnivariateConstRegion& x) {
		_a = x._a;
		_b = x._b;
		_helper = x._helper;
		_log_w_max = x._log_w_max;
		_log_w_min = x._log_w_min;
		_log_prob = x._log_prob;
		return *this;
	}
};

UnivariateConstRegion::UnivariateConstRegion(double a, double b, const UnivariateHelper<double>& helper)
: _a(a), _b(b), _helper(&helper)
{
	if (a > b) {
		Rcpp::stop("a > b");
	}

	_log_w_max = optimize(true);
	_log_w_min = optimize(false);

	// Compute g.p(b) - g.p(a) on the log scale
	_log_prob = log_sub2_exp(_helper->p(_b, true), _helper->p(_a, true));
}

double UnivariateConstRegion::d_base(const double& x, bool log) const
{
	return _helper->d(x, log);
}

std::vector<double> UnivariateConstRegion::r(unsigned int n) const
{
	// Generate a draw from $g_j$; i.e., the density $g$ truncated to this region.
	// Compute g$q((pb - pa) * u + pa) on the log scale
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	double log_pa = _helper->p(_a, true);
	const Rcpp::NumericVector& log_p = log_add2_exp(_log_prob + log(u), Rcpp::rep(log_pa, n));

	std::vector<double> out;
	for (unsigned int i = 0; i < n; i++) {
		out.push_back(_helper->q(log_p(i), true));
	}
	return out;
}

double UnivariateConstRegion::d(const double& x, bool log) const
{
	double out;
	if (!s(x)) {
		out = R_NegInf;
	} else {
		out = _helper->d(x, true) - log_sub2_exp(_helper->p(_b, true), _helper->p(_a, true));
	}
	return log ? out : exp(out);
}

bool UnivariateConstRegion::s(const double& x) const
{
	return (_a < x & x <= _b) && _helper->s(x);
}

double UnivariateConstRegion::w(const double& x, bool log) const
{
	return _helper->w(x, log);
}

double UnivariateConstRegion::w_major(const double& x, bool log) const
{
	double out = _helper->s(x) ? _log_w_max : R_NegInf;
	return log ? out : exp(out);
}

double UnivariateConstRegion::midpoint() const
{
	double out;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-Inf, Inf). Make a split at zero.
		out = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -Inf. Split based on right endpoint.
		out = _b - abs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is Inf. Split based on left endpoint.
		out = _a + std::fabs(_a) + 1;
	} else {
		out = (_a + _b) / 2;
	}

	return out;
}

std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate() const
{
	return bifurcate(midpoint());
}

std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate(const double& x) const
{
	UnivariateConstRegion r1(_a, x, *_helper);
	UnivariateConstRegion r2(x, _b, *_helper);
	return std::make_pair(r1, r2);
}

UnivariateConstRegion UnivariateConstRegion::singleton(const double& x) const
{
	return UnivariateConstRegion(x, x, *_helper);
}

bool UnivariateConstRegion::is_bifurcatable() const
{
	return true;
}

double UnivariateConstRegion::get_xi_upper(bool log) const
{
	double log_xi_upper = _log_w_max + _log_prob;
	return log ? log_xi_upper : exp(log_xi_upper);
}

double UnivariateConstRegion::get_xi_lower(bool log) const
{
	double log_xi_lower = _log_w_min + _log_prob;
	return log ? log_xi_lower : exp(log_xi_lower);
}

std::string UnivariateConstRegion::description() const
{
	char buf[32];
	sprintf(buf, "(%g, %g]", _a, _b);
	return buf;
}

void UnivariateConstRegion::print() const
{
	printf("Region<double> (%g, %g]\n", _a, _b);
}

double UnivariateConstRegion::optimize(bool maximize, bool log) const
{
	Rcpp::NumericVector log_w_endpoints = Rcpp::NumericVector::create(
		_helper->w(_a, true),
		_helper->w(_b, true)
	);
	log_w_endpoints = log_w_endpoints[!Rcpp::is_na(log_w_endpoints)];
	bool endpoint_pos_inf = Rcpp::is_true(Rcpp::any(Rcpp::is_infinite(log_w_endpoints) & log_w_endpoints > 0));
	bool endpoint_neg_inf = Rcpp::is_true(Rcpp::any(Rcpp::is_infinite(log_w_endpoints) & log_w_endpoints < 0));

	double out;

	if (maximize && endpoint_pos_inf) {
		out = R_PosInf;
	} else if (!maximize && endpoint_neg_inf) {
		out = R_NegInf;
	} else {

		RcppFunctionalUtilities::NelderMeadControl control;
		control.maxit = 100000;
		control.fnscale = maximize ? -1.0 : 1.0;

	    const RcppFunctionalUtilities::mv_function& f =
    	[&](const Rcpp::NumericVector& x) {
			double x_tx;

			// Transform to the interval (a,b]
			if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
				x_tx = x(0);
			} else if (std::isinf(_a) && _a < 0) {
				x_tx = _b*R::plogis(x(0), 0, 1, true, false);
			} else if (std::isinf(_b) && _b > 0) {
				x_tx = std::exp(x(0)) + _a;
			} else {
				x_tx = (_b - _a) * R::plogis(x(0), 0, 1, true, false) + _a;
			}

			// Call the weight function
			return _helper->w(x_tx, true);
		};

   		const Rcpp::NumericVector& init = Rcpp::NumericVector::create(0);
		const RcppFunctionalUtilities::NelderMeadResult& nm_out =
			RcppFunctionalUtilities::nelder_mead(init, f, control);

		if (nm_out.fail) {
			Rcpp::warning("Nelder-Mead: convergence status was ", nm_out.fail);
		}

		// In case the function is strictly increasing or decreasing, check the
		// objective value at the endpoints.
		if (maximize) {
			double max_lwe = Rcpp::max(log_w_endpoints);
			out = -std::max({nm_out.value, max_lwe});
		} else {
			double min_lwe = Rcpp::min(log_w_endpoints);
			out = std::min({nm_out.value, min_lwe});
		}

	}

	return log ? out : exp(out);
}

}

#endif
