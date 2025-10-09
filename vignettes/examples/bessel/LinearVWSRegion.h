#ifndef LINEAR_VWS_REGION_H
#define LINEAR_VWS_REGION_H

#include <Rcpp.h>
#include "vws.h"
#include "fntl.h"
#include "functions.h"

class LinearVWSRegion : public vws::Region<double>
{
private:
	double _a;
	double _b;
	double _lambda;
	double _nu;
	double _beta0_min;
	double _beta1_min;
	double _beta0_max;
	double _beta1_max;

public:
	LinearVWSRegion(double a, double b, double lambda, double nu);
	LinearVWSRegion(double a, double lambda, double nu);
	void init();

	double midpoint() const;
	double d_base(const double& x, bool log = false) const;
	double w(const double& x, bool log = true) const;

	std::vector<double> r(unsigned int n) const;
	double d(const double& x, bool log = false) const;

	double get_lower() const { return _a; }
	double get_upper() const { return _b; }
	double get_lambda() const { return _lambda; }
	double get_nu() const { return _nu; }
	double get_beta0_min() const { return _beta0_min; }
	double get_beta1_min() const { return _beta1_min; }
	double get_beta0_max() const { return _beta0_max; }
	double get_beta1_max() const { return _beta1_max; }

	bool s(const double& x) const {
		return _a < x && x <= _b;
	}

	bool is_bifurcatable() const {
		// Return true if the distance between a and b allows for two or more
		// integers.
		// return std::ceil(_a) <= std::floor(_b);
		return _b - _a > 1;
	}

	double w_major(const double& x, bool log = true) const;
	double xi_upper(bool log = true) const;
	double xi_lower(bool log = true) const;

	std::pair<LinearVWSRegion,LinearVWSRegion> bifurcate() const {
		double x = midpoint();
		return bifurcate(x);
	}

	std::pair<LinearVWSRegion,LinearVWSRegion> bifurcate(const double& x) const;

	std::string description() const {
		char buf[32];
		sprintf(buf, "(%g, %g]", _a, _b);
		return buf;
	}

	LinearVWSRegion singleton(const double& x) const {
		return LinearVWSRegion(x, _lambda, _nu);
	}

	void print() const {
		printf("Linear VWS Region (%g, %g]\n", _a, _b);
	}

	bool operator<(const LinearVWSRegion& x) const {
		return _a < x._a;
	}

	bool operator==(const LinearVWSRegion& x) const {
		return _a == x._a && _b == x._b;
	}

	const LinearVWSRegion& operator=(const LinearVWSRegion& x);
};

inline double LinearVWSRegion::w(const double& x, bool log) const
{
	double out = -std::lgamma(x + _nu + 1);
	return log ? out : std::exp(out);
}

inline double LinearVWSRegion::d_base(const double& x, bool log) const
{
	double lambda2 = _lambda * _lambda;
	return R::dpois(x, lambda2 / 4, log);
}

inline double LinearVWSRegion::d(const double& x, bool log) const
{
	double lambda2 = _lambda * _lambda;
	double mean = lambda2 * std::exp(_beta1_max) / 4;

	const fntl::density& f = [&](double x, bool log) {
		return R::dpois(x, mean, log);
	};

	const fntl::cdf& F = [&](double x, bool lower, bool log) {
		return R::ppois(x, mean, lower, log);
	};

	return fntl::d_trunc(x, _a, _b, f, F, log);
}

inline std::vector<double> LinearVWSRegion::r(unsigned int n) const
{
	double lambda2 = _lambda * _lambda;
	double mean = lambda2 * std::exp(_beta1_max) / 4;

	const fntl::cdf& F = [&](double x, bool lower, bool log) {
		return R::ppois(x, mean, lower, log);
	};

	const fntl::quantile& Finv = [&](double x, bool lower, bool log) {
		return R::qpois(x, mean, lower, log);
	};

	const Rcpp::NumericVector& a = Rcpp::rep(_a, n);
	const Rcpp::NumericVector& b = Rcpp::rep(_b, n);
	const auto& out = fntl::r_trunc(n, a, b, F, Finv);
	return Rcpp::as<std::vector<double>>(out);
}

inline double LinearVWSRegion::w_major(const double& x, bool log) const
{
	double out = s(x) ? _beta0_max + _beta1_max * x : R_NegInf;
	return log ? out : std::exp(out);
}

inline double LinearVWSRegion::midpoint() const
{
	double out;

	if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
		// In this case, we have an interval (-Inf, Inf). Make a split at zero.
		out = 0;
	} else if (std::isinf(_a) && _a < 0) {
		// Left endpoint is -Inf. Split based on right endpoint.
		out = _b - std::fabs(_b) - 1;
	} else if (std::isinf(_b) && _b > 0) {
		// Right endpoint is Inf. Split based on left endpoint.
		out = _a + std::fabs(_a) + 1;
	} else {
		out = (_a + _b) / 2;
	}

	return out;
}

inline LinearVWSRegion::LinearVWSRegion(double a, double b, double lambda, double nu)
: _a(a), _b(b), _lambda(lambda), _nu(nu),
	_beta0_min(), _beta1_min(), _beta0_max(), _beta1_max()
{
	init();
}

inline LinearVWSRegion::LinearVWSRegion(double a, double lambda, double nu)
: _a(a), _b(a), _lambda(lambda), _nu(nu),
	_beta0_min(), _beta1_min(), _beta0_max(), _beta1_max()
{
	// Special handling for singleton sets.
	_beta0_min = a;
	_beta1_min = 0;
	_beta0_max = a;
	_beta1_max = 0;
}

inline void LinearVWSRegion::init()
{
	// Rprintf("init checkpoint 0\n");

	if (_a >= _b) {
		Rcpp::stop("a >= b: %g >= %g", _a, _b);
	}

	if (_a == _b) {
		// Special handling for singleton sets, since truncated MGF is not
		// defined for these.
		_beta0_min = NAN;
		_beta1_min = NAN;
		_beta0_max = NAN;
		_beta1_max = NAN;
		return;
	}

	// First derivative of log-weight function
	const std::function<double(double)>& d_log_w = [&](double x) -> double
	{
		return (x + _nu + 1 > 0) ? -R::digamma(x + _nu + 1) : R_PosInf;
	};

	// Rprintf("init checkpoint 1\n");

	const std::function<double(double)>& obj = [&](double x) -> double
	{
		double dx = d_log_w(x);
		double lmgf = mgf_truncpois(dx, _lambda, _a, _b, true);
		// Rprintf("x = %g, z = %g, _lambda = %g, dx = %g, w(x, true) = %g, lmgf = %g\n",
		// 	x, _z, _lambda, dx, w(x, true), lmgf);
		return w(x, true) - x * dx + lmgf;
	};

	// Rprintf("init checkpoint 2\n");

	// log(w(x)) is concave. Any tangent line will lie above the curve.
	// Rprintf("init checkpoint 3.1.1\n");

	// For the majorizer, solve a minimization problem
	double init = midpoint();
	const auto& opt_out = vws::optimize_hybrid(obj, init, _a, _b, false);
	double c_star = opt_out.par;
	_beta0_max = w(c_star, true) - c_star*d_log_w(c_star);
	_beta1_max = d_log_w(c_star);

	// Rprintf("init checkpoint 3.1.2\n");

	// For the minorizer
	_beta1_min = (w(_b, true) - w(_a, true)) / (_b - _a);
	_beta0_min = w(_a, true) - _a*_beta1_min;

	// Rprintf("init checkpoint 3.1.3\n");

	// Rprintf("a = %g, b = %g\n", _a, _b);
	// Rprintf("beta0_min = %g, beta1_min = %g, beta0_max = %g, beta1_max = %g\n",
	// 	_beta0_min, _beta1_min, _beta0_max, _beta1_max);
}

inline double LinearVWSRegion::xi_upper(bool log) const
{
	// Compute the probability using both lower or upper tail CDF. If one is
	// unstable, it will return a -Inf and the other will be used.
	double lo = std::max(_a, 0.0);
	double hi = std::max(_b, 0.0);
	double lambda2 = _lambda * _lambda;
	double mean = lambda2 * std::exp(_beta1_max) / 4;

	double lpa = R::ppois(lo, mean, true, true);
	double lpb = R::ppois(hi, mean, true, true);
	double clpa = R::ppois(lo, mean, false, true);
	double clpb = R::ppois(hi, mean, false, true);
	double lp = vws::log_sub2_exp(lpb, lpa);
	double clp = vws::log_sub2_exp(clpa, clpb);
	double lm = std::max(lp, clp);

	// double lg1 = incgamma(std::floor(_b) + 1, mean, false, true);
	// double lg2 = incgamma(std::floor(_b) + 1, 0, false, true);
	// double lg3 = incgamma(std::ceil(_a), mean, false, true);
	// double lg4 = incgamma(std::ceil(_a), 0, false, true);
	// double lm = vws::log_sub2_exp(lg1 - lg2, lg3 - lg4);

	// Rprintf("xi_upper: a = %g, b = %g\n", _a, _b);
	// Rprintf("xi_upper: lg1 = %g, lg2 = %g, lg3 = %g, lg4 = %g, beta1_max = %g, lm = %g\n",
	// 	lg1, lg2, lg3, lg4, _beta1_max, lm);

	// double out = _beta0_max + lambda2 / 4 * std::expm1(_beta1_max) + lm;
	double out = _beta0_max - lambda2 / 4 + mean + lm;
	return log ? out : exp(out);
}

inline double LinearVWSRegion::xi_lower(bool log) const
{
	// Compute the probability using both lower or upper tail CDF. If one is
	// unstable, it will return a -Inf and the other will be used.
	double lo = std::max(_a, 0.0);
	double hi = std::max(_b, 0.0);
	double lambda2 = _lambda * _lambda;
	double mean = lambda2 * std::exp(_beta1_min) / 4;

	double lpa = R::ppois(lo, mean, true, true);
	double lpb = R::ppois(hi, mean, true, true);
	double clpa = R::ppois(lo, mean, false, true);
	double clpb = R::ppois(hi, mean, false, true);
	double lp = vws::log_sub2_exp(lpb, lpa);
	double clp = vws::log_sub2_exp(clpa, clpb);
	double lm = std::max(lp, clp);

	// double lg1 = incgamma(std::floor(_b) + 1, mean, false, true);
	// double lg2 = incgamma(std::floor(_b) + 1, 0, false, true);
	// double lg3 = incgamma(std::ceil(_a), mean, false, true);
	// double lg4 = incgamma(std::ceil(_a), 0, false, true);
	// double lm = vws::log_sub2_exp(lg1 - lg2, lg3 - lg4);

	// double out = _beta0_min + lambda2 / 4 * std::expm1(_beta1_min) + lm;
	double out = _beta0_min - lambda2 / 4 + mean + lm;

	// Ensure that xi_lower is <= xi_upper. Even if they are both coded
	// correctly, this can happen numerically. If it does happen, just take
	// lower to be equal to upper.
	double log_xi_upper = xi_upper(true);
	if (log_xi_upper < out) {
		Rprintf("LinearVWSRegion: log_xi_lower (%g) <- log_xi_upper (%g)\n",
			out, log_xi_upper);
		out = log_xi_upper;
	}

	return log ? out : exp(out);
}

inline std::pair<LinearVWSRegion,LinearVWSRegion> LinearVWSRegion::bifurcate(const double& x) const
{
	LinearVWSRegion r1(_a, x, _lambda, _nu);
	LinearVWSRegion r2(x, _b, _lambda, _nu);
	return std::make_pair(r1, r2);
}

inline const LinearVWSRegion& LinearVWSRegion::operator=(const LinearVWSRegion& x)
{
	_a = x._a;
	_b = x._b;
	_lambda = x._lambda;
	_nu = x._nu;
	_beta0_min = x._beta0_min;
	_beta1_min = x._beta1_min;
	_beta0_max = x._beta0_max;
	_beta1_max = x._beta1_max;
	return *this;
}

#endif
