// [[Rcpp::depends(vws)]]
#include "vws.h"

class MyHelper : public vws::UnivariateHelper<double>
{
public:
	MyHelper(double mu, double sigma2, double z, double lambda2)
		: _mu(mu), _sigma2(sigma2), _z(z), _lambda2(lambda2)
	{
	}

	double d(double x, bool log = false) const {
		return R::dnorm(x, _z, std::sqrt(_lambda2), log);
	}
	double p(double q, bool lower = true, bool log = false) const {
		return R::pnorm(q, _z, std::sqrt(_lambda2), lower, log);
	}
	double q(double p, bool lower = true, bool log = false) const {
		return R::qnorm(p, _z, std::sqrt(_lambda2), lower, log);
	}
	double s(double x) const {
		return 1.0;
	}
	double w(double x, bool log = false) const {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - _mu, 2.0) / (2*_sigma2);
		}
		return log ? out : exp(out);
	}
	const MyHelper& operator=(const MyHelper& x) {
		_mu = x._mu;
		_sigma2 = x._sigma2;
		_z = x._z;
		_lambda2 = x._lambda2;
		return *this;
	}

private:
	double _mu;
	double _sigma2;
	double _z;
	double _lambda2;
};

class CustomConstRegion : public vws::UnivariateConstRegion
{
private:
	double _mu;
	double _sigma2;

public:
	CustomConstRegion(double a, double b, double mu, double sigma2, const vws::UnivariateHelper<double>& helper)
	: vws::UnivariateConstRegion(a, b, helper), _mu(mu), _sigma2(sigma2)
	{
	}

	double optimize(bool maximize = true, bool log = true) const {
		double y_star = exp(_mu - _sigma2);
		double out;

		if (maximize) {
			if (y_star > _b) {
				out = _helper->w(_b, true);
			} else if (y_star < _a) {
				out = _helper->w(_a, true);
			} else {
				out = _helper->w(y_star, true);
			}
		} else {
			out = std::min(_helper->w(_a, true), _helper->w(_b, true));
		}

		return log ? out : exp(out);
	}

	// TBD: Maybe change this into a function that just finds the midpoint?
	std::pair<CustomConstRegion,CustomConstRegion> bifurcate() const {
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

		return bifurcate(x);
	}

	std::pair<CustomConstRegion,CustomConstRegion> bifurcate(const double& x) const {
		CustomConstRegion r1(_a, x, _mu, _sigma2, *_helper);
		CustomConstRegion r2(x, _b, _mu, _sigma2, *_helper);
		return std::make_pair(r1, r2);
	}

	CustomConstRegion singleton(const double& x) const {
		return CustomConstRegion(x, x, _mu, _sigma2, *_helper);
	}

	bool operator<(const CustomConstRegion& x) const {
		return UnivariateConstRegion::operator<(x);
	}

	bool operator==(const CustomConstRegion& x) const {
		return UnivariateConstRegion::operator==(x);
	}

	const CustomConstRegion& operator=(const CustomConstRegion& x) {
		UnivariateConstRegion::operator=(x);
		_mu = x._mu;
		_sigma2 = x._sigma2;
		return *this;
	}
};

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2, unsigned int N = 10, unsigned int max_rejects = 10000,
	unsigned int report_period = 1000)
{
	MyHelper helper(mu, sigma2, z, lambda2);
	CustomConstRegion supp(0.0, R_PosInf, mu, sigma2, helper);

	const std::vector<CustomConstRegion>& regions = { supp };

	vws::FMMProposal<double, CustomConstRegion> h(regions);

	h.adapt(N - 1);
	h.print(5);

	vws::RejectionControl control(max_rejects, report_period, vws::MaxRejectsAction::stop);
	const std::pair<std::vector<double>, std::vector<unsigned int>>& out =
		vws::rejection(h, n, control);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejects") = out.second
	);
}
