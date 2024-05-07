#ifndef UNIVARIATE_CONST_REGION_H
#define UNIVARIATE_CONST_REGION_H

#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <memory>
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
private:
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
	//' Majorized weight function \eqn{\overline{w}_j} for this region.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double w_major(const double& x, bool log = true) const;

	UnivariateConstRegion bifurcate_first() const;
	UnivariateConstRegion bifurcate_second() const;

	//' @description
	//' Bifurcate this region into two regions. Use \code{x} as the bifurcation
	//' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
	//' @param x An optional bifurcation point.
	UnivariateConstRegion bifurcate_first(const double& x) const;
	UnivariateConstRegion bifurcate_second(const double& x) const;

	UnivariateConstRegion singleton(const double& x) const;

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

UnivariateConstRegion UnivariateConstRegion::bifurcate_first() const
{
	double x;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-Inf, Inf). Make a split at zero.
		x = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -Inf. Split based on right endpoint.
		x = _b - abs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is Inf. Split based on left endpoint.
		x = _a + std::fabs(_a) + 1;
	} else {
		x = (_a + _b) / 2;
	}

	return bifurcate_first(x);
}

UnivariateConstRegion UnivariateConstRegion::bifurcate_second() const
{
	double x;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-Inf, Inf). Make a split at zero.
		x = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -Inf. Split based on right endpoint.
		x = _b - std::fabs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is Inf. Split based on left endpoint.
		x = _a + fabs(_a) + 1;
	} else {
		x = (_a + _b) / 2;
	}

	return bifurcate_second(x);
}

UnivariateConstRegion UnivariateConstRegion::bifurcate_first(const double& x) const
{
	return UnivariateConstRegion(_a, x, *_helper);
}

UnivariateConstRegion UnivariateConstRegion::bifurcate_second(const double& x) const
{
	return UnivariateConstRegion(x, _b, *_helper);
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

struct ext_data
{
	double a;
	double b;
	double fnscale;
	const UnivariateHelper<double>* helper;
};

double f_opt(int n, double *par, void *ex)
{
	ext_data* exd = (ext_data*) ex;

	double a = exd->a;
	double b = exd->b;
	double fnscale = exd->fnscale;
	const UnivariateHelper<double>* helper = exd->helper;
	double x = par[0];
	double x_tx;

	// Transform to the interval (a,b]
	if (std::isinf(a) && std::isinf(b) && a < 0 && b > 0) {
		x_tx = x;
	} else if (std::isinf(a) && a < 0) {
		x_tx = b*R::plogis(x, 0, 1, true, false);
	} else if (std::isinf(b) && b > 0) {
		x_tx = std::exp(x) + a;
	} else {
		x_tx = (b-a) * R::plogis(x, 0, 1, true, false) + a;
	}

	// Call the weight function
	// Rprintf("Trace: f_opt, x = %g, x_tx = %g, w(x_tx, true) = %g\n", x, x_tx, helper->w(x_tx, true));
	return fnscale * helper->w(x_tx, true);
}

double UnivariateConstRegion::optimize(bool maximize, bool log) const
{
	// control = list(maxit = 100000, warn.1d.NelderMead = false);

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

		// Rprintf("Trace: About to attempt Nelder-Mead call\n");

		// Nelder-Mead algorithm from R.
		// See https://cran.r-project.org/doc/manuals/R-exts.html#Optimization
		// and https://stackoverflow.com/questions/12765304/calling-r-function-optim-from-c
		int fail;
		int fncount;
		double f_val;
		double xin = 0;
		double par;
		double mach_eps = sqrt(std::numeric_limits<double>::epsilon());
		double fnscale = maximize ? -1.0 : 1.0;

		struct ext_data ex = { _a, _b, fnscale, _helper };

		// Rprintf("Trace: a = %f, b = %f, fnscale = %f\n", ex.a, ex.b, ex.fnscale);

		nmmin(
			1L,        // In:  int n [number of parameters]
			&xin,      // In:  double *xin [initial value]
			&par,      // Out: double *x [point at which optimum is found]
			&f_val,    // Out: double *Fmin [objective value at which optimum is found]
			f_opt,     // In:  optimfn fn [objective function]
			&fail,     // Out: int *fail [true if the function failed]
			R_NegInf,  // In:  double abstol [absolute tolerance]
			mach_eps,  // In:  double intol [user-initialized conversion tolerance]
			&ex,       // In:  void *ex [external data to pass to the objective function]
			1.0,       // In:  double alpha [reflection factor]
			0.5,       // In:  double beta [contraction and reduction factor]
			2.0,       // In:  double gamma [extension factor]
			0,         // In:  int trace [if positive, print progress info]
			&fncount,  // Out: int *fncount [number of times the objective function was called]
			100000     // In:  int maxit [maximum number of iterations]
		);

		// Rprintf("Trace: Finished Nelder-Mead call\n");

		if (fail) {
			Rcpp::warning("opt_out: convergence status was ", fail);
		}

		// In case the function is strictly increasing or decreasing, check the
		// objective value at the endpoints.
		if (maximize) {
			double x1 = std::floor(par);
			double x2 = std::ceil(par);
			double f1 = f_opt(1L, &x1, &ex);
			double f2 = f_opt(1L, &x2, &ex);
			// Rcpp::print(log_w_endpoints);
			// Rprintf("Trace: about to call std::max\n");

			double max_lwe = Rcpp::max(log_w_endpoints);
			out = -std::max({f_val, f1, f2, max_lwe});
			// Rcpp::NumericVector el = concat
			// out = std::max(f1, std::max(f2, Rcpp::max(log_w_endpoints)));
			// Rprintf("Trace: Result of optimize is f_val = %g, f1 = %g, f2 = %g, max_lwe = %g, out = %g\n",
			// 	f_val, f1, f2, max_lwe, out);
		} else {
			double x1 = std::floor(par);
			double x2 = std::ceil(par);
			double f1 = f_opt(1L, &x1, &ex);
			double f2 = f_opt(1L, &x2, &ex);
			// Rprintf("Trace: about to call std::min\n");
			double min_lwe = Rcpp::min(log_w_endpoints);
			// out = std::min(f1, f2, Rcpp::max(log_w_endpoints));
			out = std::min({f_val, f1, f2, min_lwe});
			// Rprintf("Trace: Result of optimize is f_val = %g, f1 = %g, f2 = %g, min_lwe = %g, out = %g\n",
			// 	f_val, f1, f2, min_lwe, out);
		}
	}

	return log ? out : exp(out);
}

}

#endif
