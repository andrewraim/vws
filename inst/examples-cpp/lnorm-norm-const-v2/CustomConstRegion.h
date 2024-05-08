// [[Rcpp::depends(vws)]]
#ifndef CUSTOM_CONST_REGION_H
#define CUSTOM_CONST_REGION_H

#include "vws.h"

class CustomConstRegion : public vws::UnivariateConstRegion
{
private:
	double _mu;
	double _sigma2;

public:
	CustomConstRegion(double a, double b, double mu, double sigma2,
		const vws::UnivariateHelper<double>& helper)
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

#endif
