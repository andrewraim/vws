#ifndef UNIVARIATE_CONST_REGION_FUNCTIONAL_H
#define UNIVARIATE_CONST_REGION_FUNCTIONAL_H

#include <Rcpp.h>
#include "UnivariateHelper.h"
#include "nelder-mead.h"

namespace vws {

class UnivariateConstRegionFunctional : public NelderMeadFunctional
{
protected:
	double _a;
	double _b;
	const UnivariateHelper<double>& _helper;

public:
	UnivariateConstRegionFunctional(double a, double b,
		const UnivariateHelper<double>& helper, double fnscale)
	: NelderMeadFunctional(fnscale), _helper(helper), _a(a), _b(b)
	{
	}

	virtual double operator()(const Rcpp::NumericVector& x) const {
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
		return _helper.w(x_tx, true);
	}
};

}

#endif
