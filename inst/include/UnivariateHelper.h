#ifndef UNIVARIATE_HELPER_H
#define UNIVARIATE_HELPER_H

#include "fntl.h"
#include "typedefs.h"

namespace vws {

class UnivariateHelper
{
public:
	UnivariateHelper(const fntl::density& d,
		const fntl::cdf& p, const fntl::quantile& q, const supp& s)
	: _d(d), _p(p), _q(q), _s(s)
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
		return _s(x);
	}
	const UnivariateHelper& operator=(const UnivariateHelper& x) {
		_d = x._d;
		_p = x._p;
		_q = x._q;
		_s = x._s;
		return *this;
	}

private:
	fntl::density _d;
	fntl::cdf _p;
	fntl::quantile _q;
	supp _s;
};

}

#endif
