// [[Rcpp::depends(vws)]]
#ifndef CUSTOM_LINEAR_REGION_H
#define CUSTOM_LINEAR_REGION_H

#include "vws.h"
#include "normal_truncated.h"

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
	CustomLinearRegion::CustomLinearRegion(double a, double b, double mu,
		double sigma2, double z, double lambda2);

	double mgf(double s, bool log = false, double tol = 1e-6);

	// First derivative of log w(x)
	double d_log_w(double x) {
		return -1/x * (1 + (std::log(x) - _mu) / _sigma2);
	}

	double d_base(const double& x, bool log = false) const
	{
		return dnorm(x, _z, std::sqrt(_lambda2), log);
	}

	double w(const double& x, bool log = false) const {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - _mu, 2.0) / (2*_sigma2);
		}
		return log ? out : exp(out);
	}

	std::vector<double> r(unsigned int n) const {
		double mean = _z + _beta1_max * _lambda2
		double sd = std::sqrt(_lambda2)
		r_norm_trunc(n, mean, sd, _a, _b);
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

	double xi_upper(bool log = true, double tol = 1e-6) const;
	double xi_lower(bool log = true, double tol = 1e-6) const;

	std::pair<CustomLinearRegion,CustomLinearRegion> bifurcate() const {
		if (std::isinf(_a) && std::isinf(_b) && _a < 0 && _b > 0) {
			// Here we have an interval (-Inf, Inf). Make a split at zero.
			x = 0;
		} else if (std::isinf(_a) && _a < 0) {
			// Left endpoint is -Inf. Split based on right endpoint.
			x = b - abs(b) - 1;
		} else if (std::isinf(_b) && _b > 0) {
			// Right endpoint is Inf. Split based on left endpoint.
			x = a + abs(_a) + 1;
		} else {
			x = (_a + _b) / 2;
		}

		return bifurcate(x);
	}

	std::pair<CustomLinearRegion,CustomLinearRegion> bifurcate(const double& x) const {
		CustomLinearRegion r1(_a, x, _mu, _sigma2, _z, _lambda2);
		CustomLinearRegion r2(x, _b, _mu, _sigma2, _z, _lambda2);
		return std::make_pair(r1, r2);
	}

	std::string description() const {
		charbuf[32];
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
		_a = a;
		_b = b;
		_mu = mu;
		_sigma2 = sigma2;
		_z = z;
		_lambda2 = lambda2;
		_beta0_min = beta0_min;
		_beta1_min = beta1_min;
		_beta0_max = beta0_max;
		_beta1_max = beta1_max;
		return *this;
	}
};

// MGF the truncated and reweighted g
double CustomLinearRegion::mgf(double s, bool log = false, double tol)
{
	// If we are truncating way into the upper tail of a distribution, working
	// with the complement of the CDF helps to retain precision. Otherwise,
	// work with the CDF function.

	double lp_num_a = R::pnorm(_a, _z + s*_lambda2, std::sqrt(_lambda2), true, true);
	double lp_num_b = R::pnorm(_b, _z + s*_lambda2, std::sqrt(_lambda2), true, true);
	double clp_num_a = R::pnorm(_a, _z + s*_lambda2, std::sqrt(_lambda2), false, true);
	double clp_num_b = R::pnorm(_b, _z + s*_lambda2, std::sqrt(_lambda2), false, true);
	double lp_num = lp_num_a > log1p(-tol) ? vws::log_sub2_exp(clp_num_a, clp_num_b) : vws::log_sub2_exp(lp_num_b, lp_num_a);

	double lp_den_a = R::pnorm(_a, _z, std::sqrt(_lambda2), true, true);
	double lp_den_b = R::pnorm(_b, _z, std::sqrt(_lambda2), true, true);
	double clp_den_a = R::pnorm(_a, _z, std::sqrt(_lambda2), false, true);
	double clp_den_b = R::pnorm(_b, _z, std::sqrt(_lambda2), false, true);
	double lp_den = lp_den_a > log1p(-tol) ? vws::log_sub2_exp(clp_den_a, clp_den_b) : vws::log_sub2_exp(lp_den_b, lp_den_a);

	double out = lp_num - lp_den + s*_z + pow(s, 2.0) * _lambda2 / 2;
	return log ? out : exp(out);
}

CustomLinearRegion::CustomLinearRegion(double a, double b, double mu,
	double sigma2, double z, double lambda2)
: _a(a), _b(b), _mu(mu), _sigma2(sigma2), _z(z), _lambda2(lambda2),
  _beta0_min(), _beta1_min(), _beta0_max(), _beta1_max()
{
	if (a > b) {
		Rcpp::stop("a > b");
	}

	obj_line = function(x) {
		double gr = d_log_w(x);
		return w(x, true) - x * gr + mgf(gr, true);
	}

	bool l_concave = log(a) < mu - sigma2 + 1;
	bool r_convex = log(b) > mu - sigma2 + 1;
	if (l_concave && r_convex) {
		Rcpp::stop("Partition your region so that %g is not in the interior\n", exp(mu - sigma2 + 1));
	}

	is_concave = l_concave;
	is_convex = r_convex;

	if (is_concave) {
		// log w(x) is concave

		// For the minorizer
		optim_out = optimize(f = obj_line, interval = c(a, b), maximum = FALSE);
		c_star = optim_out$minimum;
		_beta0_max = w(c_star) - c_star*d_log_w(c_star);
		_beta1_max = d_log_w(c_star);

		// For the majorizer
		A = matrix(c(1,1,a,b), 2, 2);
		c = c(w(a), w(b));
		x = solve(A, c);
		_beta0_min = x[0];
		_beta1_min = x[1];
	} else if (is_convex) {
		// log w(x) is convex

		// For the minorizer
		optim_out = optimize(f = obj_line, interval = c(a, b), maximum = FALSE);
		c_star = optim_out$minimum;
		_beta0_min = w(c_star) - c_star*d_log_w(c_star);
		_beta1_min = d_log_w(c_star);

		// For the majorizer
		A = matrix(c(1,1,a,b), 2, 2);
		c = c(w(a), w(b));
		x = solve(A, c);
		_beta0_max = x[0];
		_beta1_max = x[1];
	} else {
		// log w(x) is constant
		_beta0_min = 0;
		_beta1_min = 0;
		_beta0_max = 0;
		_beta1_max = 0;
	}
}

double CustomLinearRegion::xi_upper(bool log = true, double tol) const
{
	double out;

	double lp_a = R::pnorm(_a, _z + _beta1_max*_lambda2, sqrt(_lambda2), true, true);
	if (lp_a > std::log1p(-tol)) {
		// If we are truncating way into the upper tail of a distribution, working
		// with the complement of the CDF helps to retain precision.
		double clp_a = pnorm(_a, _z + _beta1_max*_lambda2, std::sqrt(_lambda2), false, true);
		double clp_b = pnorm(_b, _z + _beta1_max*_lambda2, std::sqrt(_lambda2), false, true);
		out = _beta0_max + _beta1_max * _z + _beta1_max^2 * _lambda2 / 2 +
			vws::log_sub2_exp(clp_a, clp_b)
	} else {
		// Otherwise, work with the CDF function.
		double lp_b = pnorm(_b, _z + _beta1_max*_lambda2, std::sqrt(_lambda2), true, true)
		out = _beta0_max + _beta1_max * _z + _beta1_max^2 * _lambda2 / 2 +
			vws::log_sub2_exp(lp_b, lp_a)
	}

	return log ? out : exp(out);
}

double xi_lower(bool log = true, double tol) const
{
	double lp_a = R::pnorm(_a, _z + _beta1_min*_lambda2, sqrt(_lambda2), true, true);
	if (lp_a > std::log1p(-tol)) {
		// If we are truncating way into the upper tail of a distribution, working
		// with the complement of the CDF helps to retain precision.
		double clp_a = pnorm(_a, _z + _beta1_min*_lambda2, std::sqrt(_lambda2), false, true);
		double clp_b = pnorm(_b, _z + _beta1_min*_lambda2, std::sqrt(_lambda2), false, true);
		double out = _beta0_min + _beta1_min * _z + _beta1_min^2 * _lambda2 / 2 +
			vws::log_sub2_exp(clp_a, clp_b);
	} else {
		// Otherwise, work with the CDF function.
		double lp_b = pnorm(_b, _z + _beta1_min*_lambda2, std::sqrt(_lambda2), true, true);
		double out = _beta0_min + _beta1_min * _z + _beta1_min^2 * _lambda2 / 2 +
			vws::log_sub2_exp(lp_b, lp_a);
	}

	// If we hopelessly run out of precision and the result is larger than
	// xi_upper, set xi_lower to xi_upper.
	out = std::min(out, xi_upper(true, tol));
	return log ? out : exp(out);
}

#endif
