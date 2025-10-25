#ifndef UNIVARIATE_HELPER_H
#define UNIVARIATE_HELPER_H

#include "fntl.h"
#include "typedefs.h"

namespace vws {

class UnivariateHelper
{
public:
	UnivariateHelper(const fntl::density& d, const fntl::cdf& p,
		const fntl::quantile& q)
	: _d(d), _p(p), _q(q)
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
	const UnivariateHelper& operator=(const UnivariateHelper& x) {
		_d = x._d;
		_p = x._p;
		_q = x._q;
		return *this;
	}

private:
	fntl::density _d;
	fntl::cdf _p;
	fntl::quantile _q;
};

}

#endif

