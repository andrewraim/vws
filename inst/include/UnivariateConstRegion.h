#ifndef UNIVARIATE_CONST_H
#define UNIVARIATE_CONST_H

#include <Rcpp.h>
#include <R_ext/Applic.h>
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
	UnivariateHelper<double> _g;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;

public:
	virtual double w(double x, bool log = false) = 0;

	//' @param a Lower limit of interval.
	//' @param b Upper limit of interval.
	//' @param w Weight function for the target distribution.
	//' @param g An object created by \code{univariate_helper}.
	UnivariateConstRegion(double a, double b, const UnivariateHelper<double>& g);

	//' @description
	//' Density function \eqn{g} for the base distribution.
	//' @param x Density argument.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double d_base(double x, bool log = false) const;

	//' @description
	//' Generate a draw from \eqn{g_j} specific to this region.
	//' @param n Number of draws to generate.
	//' @return A list of draws, with one draw per list element.
	std::vector<double> r(unsigned int n) const;

	//' @description
	//' Density of \eqn{g_j} specific to this region.
	//' @param x Density argument.
	double d(double x) const;

	//' @description
	//' Test if given \code{x} is in the support for the \eqn{g_j} specific to
	//' this region.
	//' @param x Density argument.
	bool s(double x) const;

	//' @description
	//' Majorized weight function \eqn{\overline{w}_j} for this region.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double w_major(double x, bool log = true) const;

	std::pair<UnivariateConstRegion,UnivariateConstRegion> bifurcate() const;

	//' @description
	//' Bifurcate this region into two regions. Use \code{x} as the bifurcation
	//' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
	//' @param x An optional bifurcation point.
	std::pair<UnivariateConstRegion,UnivariateConstRegion> bifurcate(double x) const;

	//' @description
	//' Return a logical value indicating whether this region is bifurcatable.
	bool is_bifurcatable() const;

	//' @description
	//' The quantity \eqn{\overline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double xi_upper(bool log = true) const;

	//' @description
	//' The quantity \eqn{\underline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double xi_lower(bool log = true) const;

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
};

UnivariateConstRegion::UnivariateConstRegion(double a, double b, const UnivariateHelper<double>& g)
: _a(a), _b(b), _g(g)
{
	if (a > b) {
		Rcpp::stop("a > b");
	}

	_log_w_max = optimize(true);
	_log_w_min = optimize(false);

	// Compute g.p(b) - g.p(a) on the log scale
	_log_prob = log_sub2_exp(g.p(_b, true), g.p(_a, true));
}


double UnivariateConstRegion::d_base(double x, bool log = false) const
{
	return _g.d(x, log);
}

std::vector<double> UnivariateConstRegion::r(unsigned int n) const
{
	// Generate a draw from $g_j$; i.e., the density $g$ truncated to this region.
	// Compute g$q((pb - pa) * u + pa) on the log scale
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	double log_pa = _g.p(_a, true);
	const Rcpp::NumericVector& log_p = log_add2_exp(_log_prob + log(u), Rcpp::rep(log_pa, n));
	return _g.q(log_p, true);
}

double UnivariateConstRegion::d(double x) const
{
	double out;
	if (!s(x)) {
		out = R_NegInf;
	} else {
		out = _g.d(x, true) - log_sub2_exp(_g.p(_b, true), _g.p(_a, true));
	}
	return log ? out : exp(out);
}

bool UnivariateConstRegion::s(double x) const
{
	return (_a < x & x <= _b) && _g.s(x);
}

double UnivariateConstRegion::w_major(double x, bool log = true) const
{
	double out = _g.s(x) ? _log_w_max : R_NegInf;
	return log ? out : exp(out);
}

std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate() const
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
		x = _a + fabs(_a) + 1;
	} else {
		x = (_a + _b) / 2;
	}

	return bifurcate(x);
}

std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate(double x) const
{
	UnivariateConstRegion s1(_a, x, _g);
	UnivariateConstRegion s2(x, _b, _g);
	std::make_pair(s1, s2);
}

bool UnivariateConstRegion::is_bifurcatable() const
{
	return true;
}

double UnivariateConstRegion::xi_upper(bool log = true) const
{
	double log_xi_upper = _log_w_max + _log_prob;
	return log ? log_xi_upper : exp(log_xi_upper);
}

double UnivariateConstRegion::xi_lower(bool log = true) const
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
	printf("Univariate Const Region (%g, %g]\n", _a, _b);
}

double UnivariateConstRegion::optimize(bool maximize = true, bool log = true) const
{
	// control = list(maxit = 100000, warn.1d.NelderMead = false);

	double f_opt(int n, double *par, void *ex)
	{
		double a = ex->a;
		double b = ex->b;
		double fnscale = ex->fnscale;
		double x = par[0];

		// Transform to the interval (a,b]
		if (Rcpp::is_infinite(a) && Rcpp::is_infinite(b) && a < 0 && b > 0) {
			x_tx = x;
		} else if (Rcpp::is_infinite(a) && a < 0) {
			x_tx = b*plogis(x);
		} else if (Rcpp::is_infinite(b) && b > 0) {
			x_tx = Rcpp::exp(x) + a;
		} else {
			x_tx = (b-a) * Rcpp::plogis(x) + a;
		}

		// Call the weight function
		fnscale * w(x_tx, true);
	}

	struct ext_data
	{
		int a;
		int b;
		double fnscale;
	};

	Rcpp::LogicalVector log_w_endpoints = Rcpp::LogicalVector::create(w(_a, true), w(_b, true));
	log_w_endpoints = log_w_endpoints[!Rcpp::is_na(log_w_endpoints)];
	bool endpoint_pos_inf = Rcpp::any(Rcpp::is_infinite(log_w_endpoints) & log_w_endpoints > 0);
	bool endpoint_neg_inf = Rcpp::any(Rcpp::is_infinite(log_w_endpoints) & log_w_endpoints < 0);

	double out;

	if (maximize && endpoint_pos_inf) {
		out = R_PosInf;
	} else if (!maximize && endpoint_neg_inf) {
		out = R_NegInf;
	} else {
		// Nelder-Mead algorithm from R.
		// See https://cran.r-project.org/doc/manuals/R-exts.html#Optimization
		// and https://stackoverflow.com/questions/12765304/calling-r-function-optim-from-c
		int fail;
		int fncount;
		double Fmin;
		double xin = 0;
		double par;
		double mach_eps = sqrt(std::numeric_limits<double>::epsilon);

		struct ext_data ex = { _a, _b, -1*maximize + 1*minimize };

		nmmin(
			1L,        // In:  int n [number of parameters]
			&xin,      // In:  double *xin [initial value]
			&par,      // Out: double *x [point at which optimum is found]
			&Fmin,     // Out: double *Fmin [objective value at which optimum is found]
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

		if (fail) {
			Rcpp::warning("opt_out: convergence status was ", fail);
		}

		if (maximize) {
			out = std::max(f_opt(std::floor(par)), f_opt(std::ceil(par)), log_w_endpoints);
		} else {
			out = std::min(f_opt(std::floor(par)), f_opt(std::ceil(par)), log_w_endpoints);
		}
	}

	return log ? out : exp(out);
}

}

#endif
