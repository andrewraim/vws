// [[Rcpp::depends(vws)]]
#ifndef CUSTOM_LINEAR_REGION_H
#define CUSTOM_LINEAR_REGION_H

#include <Rcpp.h>
#include "vws.h"
#include "normal-truncated.h"
#include "nelder-mead.h"

class CustomLinearRegion : public vws::Region<double>
{
private:
	double _a;
	double _b;
	double _mu;
	double _sigma2;
	double _z;
	double _lambda2;
	double _beta0_min;
	double _beta1_min;
	double _beta0_max;
	double _beta1_max;

public:
	CustomLinearRegion(double a, double b, double mu,
		double sigma2, double z, double lambda2);

	double mgf(double s, bool log = true) const;
	double midpoint() const;

	double d_base(const double& x, bool log = false) const
	{
		return R::dnorm(x, _z, std::sqrt(_lambda2), log);
	}

	double w(const double& x, bool log = true) const {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - _mu, 2.0) / (2*_sigma2);
		}
		return log ? out : exp(out);
	}

	std::vector<double> r(unsigned int n) const {
		double mean = _z + _beta1_max * _lambda2;
		double sd = std::sqrt(_lambda2);
		const Rcpp::NumericVector& out = r_norm_trunc(n, mean, sd, _a, _b);
		return Rcpp::as<std::vector<double>>(out);
	}

	double d(const double& x, bool log = false) const {
		double mean = _z + _beta1_max * _lambda2;
		double sd = std::sqrt(_lambda2);
		return d_norm_trunc(x, mean, sd, _a, _b, log);
	}

	bool s(const double& x) const {
		return _a < x & x <= _b;
	}

	bool is_bifurcatable() const {
		return true;
	}

	double w_major(const double& x, bool log = true) const {
		double out = R_NegInf;
		if (s(x)) {
			out = _beta0_max + _beta1_max*x;
		}
		return log ? out : exp(out);
	}

	double get_xi_upper(bool log = true) const;
	double get_xi_lower(bool log = true) const;

	std::pair<CustomLinearRegion,CustomLinearRegion> bifurcate() const {
		return bifurcate(midpoint());
	}

	std::pair<CustomLinearRegion,CustomLinearRegion> bifurcate(const double& x) const {
		// printf("About to begin bifurcate at %g\n", x);
		CustomLinearRegion r1(_a, x, _mu, _sigma2, _z, _lambda2);
		CustomLinearRegion r2(x, _b, _mu, _sigma2, _z, _lambda2);
		// printf("About to return from bifurcate at %g\n", x);
		return std::make_pair(r1, r2);
	}

	std::string description() const {
		char buf[32];
		sprintf(buf, "(%g, %g]", _a, _b);
		return buf;
	}

	CustomLinearRegion singleton(const double& x) const {
		return CustomLinearRegion(x, x, _mu, _sigma2, _z, _lambda2);
	}

	void print() const {
		printf("Custom Linear Region (%g, %g]\n", _a, _b);
	}

	bool operator<(const CustomLinearRegion& x) const {
		return _a < x._a;
	}

	bool operator==(const CustomLinearRegion& x) const {
		return _a == x._a && _b == x._b;
	}

	const CustomLinearRegion& operator=(const CustomLinearRegion& x) {
		_a = x._a;
		_b = x._b;
		_mu = x._mu;
		_sigma2 = x._sigma2;
		_z = x._z;
		_lambda2 = x._lambda2;
		_beta0_min = x._beta0_min;
		_beta1_min = x._beta1_min;
		_beta0_max = x._beta0_max;
		_beta1_max = x._beta1_max;
		return *this;
	}
};

// MGF of the truncated and reweighted g
double CustomLinearRegion::mgf(double s, bool log) const
{
	if (_a >= _b) {
		return NAN;
	}

	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision. Otherwise,
	// work with the CDF function.

	double lp_num_a = R::pnorm(_a, _z + s*_lambda2, std::sqrt(_lambda2), true, true);
	double lp_num_b = R::pnorm(_b, _z + s*_lambda2, std::sqrt(_lambda2), true, true);
	double clp_num_a = R::pnorm(_a, _z + s*_lambda2, std::sqrt(_lambda2), false, true);
	double clp_num_b = R::pnorm(_b, _z + s*_lambda2, std::sqrt(_lambda2), false, true);
	double lp_num = std::max(
		vws::log_sub2_exp(clp_num_a, clp_num_b),
		vws::log_sub2_exp(lp_num_b, lp_num_a)
	);

	double lp_den_a = R::pnorm(_a, _z, std::sqrt(_lambda2), true, true);
	double lp_den_b = R::pnorm(_b, _z, std::sqrt(_lambda2), true, true);
	double clp_den_a = R::pnorm(_a, _z, std::sqrt(_lambda2), false, true);
	double clp_den_b = R::pnorm(_b, _z, std::sqrt(_lambda2), false, true);
	double lp_den = std::max(
		vws::log_sub2_exp(clp_den_a, clp_den_b),
		vws::log_sub2_exp(lp_den_b, lp_den_a)
	);

	// printf("In mgf\n");
	// printf("lp_num: %g\n", lp_num);
	// printf("lp_den: %g\n", lp_den);
	// printf("s: %g\n", s);
	// printf("_z: %g\n", _z);
	// printf("_a: %g\n", _a);
	// printf("_b: %g\n", _b);
	// printf("_lambda2: %g\n", _lambda2);

	double out = lp_num - lp_den + s*_z + pow(s, 2.0) * _lambda2 / 2;
	return log ? out : exp(out);
}

double CustomLinearRegion::midpoint() const
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

CustomLinearRegion::CustomLinearRegion(double a, double b, double mu,
	double sigma2, double z, double lambda2)
: _a(a), _b(b), _mu(mu), _sigma2(sigma2), _z(z), _lambda2(lambda2),
  _beta0_min(), _beta1_min(), _beta0_max(), _beta1_max()
{
	if (a > b) {
		Rcpp::stop("a > b");
	}

	if (a == b) {
		// Special handling for singleton sets, since truncated MGF is not
		// defined for these.
		_beta0_min = NAN;
		_beta1_min = NAN;
		_beta0_max = NAN;
		_beta1_max = NAN;
		return;
	}

	// printf("Begin constructor for CustomLinearRegion\n");

    std::function<double(double)> d_log_w = [&](double x) {
		return -1.0/x * (1.0 + (std::log(x) - _mu) / _sigma2);
   	};

    const std::function<double(double)>& tx = [&](double x) {
		// Transform to the interval (a,b]
		if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
			return x;
		} else if (std::isinf(_a) && _a < 0) {
			return _b*R::plogis(x, 0, 1, true, false);
		} else if (std::isinf(_b) && _b > 0) {
			return std::exp(x) + _a;
		} else {
			return (_b - _a) * R::plogis(x, 0, 1, true, false) + _a;
		}
    };

	RcppFunctionalUtilities::mv_function f = [&](const Rcpp::NumericVector& x) {
		double x_tx = tx(x(0));
		double gr = d_log_w(x_tx);
		return w(x_tx, true) - x_tx * gr + mgf(gr, true);
	};

	bool l_concave = std::log(_a) < _mu - _sigma2 + 1;
	bool r_convex = std::log(_b) > _mu - _sigma2 + 1;
	if (l_concave && r_convex) {
		Rcpp::stop("Partition your region so that %g is not in the interior\n", exp(_mu - _sigma2 + 1));
	}

	bool is_concave = l_concave;
	bool is_convex = r_convex;

	if (is_concave) {
		// log w(x) is concave

		// For the minorizer
		RcppFunctionalUtilities::NelderMeadControl control;
		control.maxit = 100000;
		control.fnscale = 1.0;
		const Rcpp::NumericVector& init = Rcpp::NumericVector::create(0);
		const RcppFunctionalUtilities::NelderMeadResult& nm_out = RcppFunctionalUtilities::nelder_mead(init, f, control);
		double c_star = tx(nm_out.par(0));
		_beta0_max = w(c_star) - c_star * d_log_w(c_star);
		_beta1_max = d_log_w(c_star);

		// For the majorizer
		_beta1_min = (w(b) - w(a)) / (b - a);
		_beta0_min = w(a) - a*_beta1_min;

	} else if (is_convex) {
		// log w(x) is convex

		// For the minorizer
		RcppFunctionalUtilities::NelderMeadControl control;
		control.maxit = 100000;
		control.fnscale = 1.0;
		const Rcpp::NumericVector& init = Rcpp::NumericVector::create(0);
		const RcppFunctionalUtilities::NelderMeadResult& nm_out = RcppFunctionalUtilities::nelder_mead(init, f, control);
		double c_star = tx(nm_out.par(0));
		_beta0_min = w(c_star) - c_star*d_log_w(c_star);
		_beta1_min = d_log_w(c_star);

		// For the majorizer
		_beta1_max = (w(b) - w(a)) / (b - a);
		_beta0_max = w(a) - a*_beta1_max;
	} else {
		// log w(x) is constant
		_beta0_min = 0;
		_beta1_min = 0;
		_beta0_max = 0;
		_beta1_max = 0;
	}

	if (std::isnan(_beta1_max)) {
		printf("is_convex = %d\n", is_convex);
		printf("_a = %g\n", _a);
		printf("_b = %g\n", _b);
		printf("w(_a) = %g\n", w(_a));
		printf("w(_b) = %g\n", w(_b));
		printf("_beta0_min = %g\n", _beta0_min);
		printf("_beta1_min = %g\n", _beta1_min);
		printf("_beta0_max = %g\n", _beta0_max);
		printf("_beta1_max = %g\n", _beta1_max);
		Rcpp::stop("PAUSE!");
	}

	// printf("End constructor for CustomLinearRegion\n");
}

double CustomLinearRegion::get_xi_upper(bool log) const
{
	// printf("Begin get_xi_upper\n");
	// printf("_a = %g, _b = %g, _z = %g, _lambda2 = %g\n", _a, _b, _z, _lambda2);
	// printf("_beta0_max = %g, _beta1_max = %g\n", _beta0_max, _beta1_max);

	double lp_a = R::pnorm(_a, _z + _beta1_max*_lambda2, sqrt(_lambda2), true, true);
	double lp_b = R::pnorm(_b, _z + _beta1_max*_lambda2, std::sqrt(_lambda2), true, true);
	double lp1_diff = vws::log_sub2_exp(lp_b, lp_a);

	double clp_a = R::pnorm(_a, _z + _beta1_max*_lambda2, std::sqrt(_lambda2), false, true);
	double clp_b = R::pnorm(_b, _z + _beta1_max*_lambda2, std::sqrt(_lambda2), false, true);
	double lp2_diff = vws::log_sub2_exp(clp_a, clp_b);

	// printf("lp_a = %g, lp_b = %g, lp1_diff = %g\n", lp_a, lp_b, lp1_diff);
	// printf("clp_a = %g, clp_b = %g, lp2_diff = %g\n", clp_a, clp_b, lp2_diff);

	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision.
	double lp_diff = std::max(lp1_diff, lp2_diff);
	// printf("lp_diff = %g\n", lp_diff);
	double out = _beta0_max + _beta1_max * _z + std::pow(_beta1_max, 2.0) * _lambda2 / 2 + lp_diff;

	// printf("End get_xi_upper with (log) out = %g\n", out);
	return log ? out : exp(out);
}

double CustomLinearRegion::get_xi_lower(bool log) const
{
	double lp_a = R::pnorm(_a, _z + _beta1_min*_lambda2, sqrt(_lambda2), true, true);
	double lp_b = R::pnorm(_b, _z + _beta1_min*_lambda2, std::sqrt(_lambda2), true, true);
	double lp1_diff = vws::log_sub2_exp(lp_b, lp_a);

	double clp_a = R::pnorm(_a, _z + _beta1_min*_lambda2, std::sqrt(_lambda2), false, true);
	double clp_b = R::pnorm(_b, _z + _beta1_min*_lambda2, std::sqrt(_lambda2), false, true);
	double lp2_diff = vws::log_sub2_exp(clp_a, clp_b);

	double lp_diff = std::max(lp1_diff, lp2_diff);

	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision.
	double out = _beta0_min + _beta1_min * _z + std::pow(_beta1_min, 2.0) * _lambda2 / 2 + lp_diff;

	// If we hopelessly run out of precision and the result is larger than
	// xi_upper, set xi_lower to xi_upper.
	out = std::min(out, get_xi_upper(true));
	// printf("End get_xi_lower with (log) out = %g\n", out);
	return log ? out : exp(out);
}

#endif
