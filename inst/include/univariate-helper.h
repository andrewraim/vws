#ifndef UNIVARIATE_HELPER_H
#define UNIVARIATE_HELPER_H

#include "fntl.h"
#include "typedefs.h"

namespace vws {

class univariate_helper
{
public:
	univariate_helper(fntl::density d, fntl::cdf p, fntl::quantile q)
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
	const univariate_helper& operator=(const univariate_helper& x) {
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

