#include <Rcpp.h>
#include "Region.h"
#include "UnivariateHelper.h"

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
	UnivariateHelper _g;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;

public:
	double w(double x, bool log = false) = 0;

	//' @param a Lower limit of interval.
	//' @param b Upper limit of interval.
	//' @param w Weight function for the target distribution.
	//' @param g An object created by \code{univariate_helper}.
	UnivariateConstRegion(a, b, w, g);

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

UnivariateConstRegion(double a, double b, w, const UnivariateHelper& g)
: _a(a), _b(b), _g(g), _w(w)
{
	stopifnot(a <= b);

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
	const RcppNumericVector& u = Rcpp::runif(n);
	log_pa = _g.p(_a, true);
	log_p = log_add2_exp(_log_prob + log(u), Rcpp::rep(log_pa, n));
	return _g.q(log_p, true);
}

double UnivariateConstRegion::d(double x) const
{
	double out = s(x) ?
		_g.d(x, true) - log_sub2_exp(_g.p(_b, true), _g.p(_a, true)) :
		Rcpp::NegInf;
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

	if (is.infinite(_a) && is.infinite(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-Inf, Inf). Make a split at zero.
		x = 0;
	} else if (is.infinite(_a) && _a < 0) {
		// Left endpoint is -Inf. Split based on right endpoint.
		x = b - abs(b) - 1;
	} else if (is.infinite(_b) && _b > 0) {
		// Right endpoint is Inf. Split based on left endpoint.
		x = _a + abs(_a) + 1;
	} else {
		x = (_a + _b) / 2;
	}

	return bifurcate(x);
}

std::pair<UnivariateConstRegion,UnivariateConstRegion>
UnivariateConstRegion::bifurcate(double x) const
{
	UnivariateConstRegion(_a, x, w, _g) s1;
	UnivariateConstRegion(x, _b, w, _g) s2;
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
	sprintf("(%g, %g]", _a, _b);
}

void UnivariateConstRegion::print() const
{
	printf("Univariate Const Region (%g, %g]\n", _a, _b);
}

double UnivariateConstRegion::optimize(bool maximize = true, bool log = true) const
{
	a = _a;
	b = _b;
	w = _w;

	method = "Nelder-Mead";
	control = list(maxit = 100000, warn.1d.NelderMead = false);

	// Transformation to bounded interval, if necessary
	tx = function(x) {
		if (is.infinite(a) && is.infinite(b) && a < 0 && b > 0) {
			out = x
		} else if (is.infinite(a) && a < 0) {
			out = b*plogis(x)
		} else if (is.infinite(b) && b > 0) {
			out = exp(x) + a
		} else {
			out = (b-a) * plogis(x) + a
		}
		return(out)
	}

	f_opt = function(x) {
		w(tx(x), true);
	}

	init = 0;
	log_w_endpoints = c(w(a, log = TRUE), w(b, log = TRUE));
	log_w_endpoints = log_w_endpoints[!is.na(log_w_endpoints)];

	endpoint_pos_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints > 0);
	endpoint_neg_inf = any(is.infinite(log_w_endpoints) & log_w_endpoints < 0);

	// Nelder-Mead algorithm from R
	void nmmin(
		int n,               // Number of parameters
		double *xin,         // Initial value
		double *x,           // Point at which optimum is found
		double *Fmin,        // Value at which optimum is found
		optimfn fn,          // Function to optimize
		int *fail,           // ...
		double abstol,
		double intol,
		void *ex,
		double alpha,
		double beta,
		double gamma,
		int trace,
		int *fncount,
		int maxit
	);

	if (maximize && endpoint_pos_inf) {
		out = Inf;
	} else if (maximize) {
		control$fnscale = -1;
		opt_out = optim(init, f_opt, method = method, control = control);
		if (opt_out.convergence != 0) {
			warning("opt_out: convergence status was ", opt_out.convergence);
			browser();
		}
		out = max(f_opt(floor(opt_out.par)), f_opt(ceiling(opt_out.par)), log_w_endpoints);
	}

	if (!maximize && endpoint_neg_inf) {
		out = -Inf;
	} else if (!maximize) {
		control$fnscale = 1;
		opt_out = optim(init, f_opt, method = method, control = control);
		if (opt_out.convergence != 0) {
			warning("opt_out: convergence status was ", opt_out$convergence)
			browser();
		}
		out = min(f_opt(floor(opt_out$par)), f_opt(ceiling(opt_out$par)), log_w_endpoints);
	}

	return log ? out : exp(out);
}
