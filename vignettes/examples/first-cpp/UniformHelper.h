#ifndef UNIFORM_HELPER_H
#define UNIFORM_HELPER_H

#include "vws.h"

class UniformHelper : public vws::UnivariateHelper<double>
{
public:
    UniformHelper(double a, double b)
    : _a(a), _b(b)
    {
    }

    double d(double x, bool log = false) const {
        return R::dunif(x, _a, _b, log);
    }
    double p(double q, bool lower = true, bool log = false) const {
        return R::punif(q, _a, _b, lower, log);
    }
    double q(double p, bool lower = true, bool log = false) const {
        return R::qunif(p, _a, _b, lower, log);
    }
    bool s(double x) const {
        return _a <= x && x <= _b;
    }
    const UniformHelper& operator=(const UniformHelper& x) {
        _a = x._a;
        _b = x._b;
        return *this;
    }

private:
    double _a;
    double _b;
};

#endif
