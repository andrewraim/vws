#ifndef LAMBDA_HELPER_H
#define LAMBDA_HELPER_H

#include "fntl.h"
#include "typedefs.h"
#include "UnivariateHelper.h"

/*
* Subclass of `UnivariateHelper` for functions specified via Rcpp.
*/

namespace vws {

class LambdaHelper : public UnivariateHelper<double>
{
public:
	LambdaHelper(double a, double b, const fntl::density& d,
		const fntl::cdf& p, const fntl::quantile& q, const supp& s)
	: _a(a), _b(b), _d(d), _p(p), _q(q), _s(s)
	{
	}

	double d(double x, bool log = false) const {
		return _d(x, log);
	}
	double p(double q, bool lower = true, bool log = false) const {
		return _p(q, lower, log);
	}
	double q(double p, bool lower = true, bool log = false) const {
		return _q(p, lower, log);
	}
	bool s(double x) const {
		return _a <= x && x <= _b && _s(x);
	}
	const LambdaHelper& operator=(const LambdaHelper& x) {
		_a = x._a;
		_b = x._b;
		_d = x._d;
		_p = x._p;
		_q = x._q;
		_s = x._s;
		return *this;
	}

private:
	double _a;
	double _b;
	fntl::density _d;
	fntl::cdf _p;
	fntl::quantile _q;
	supp _s;
};

}

#endif
