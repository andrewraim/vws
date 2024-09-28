#ifndef NORMAL_HELPER_H
#define NORMAL_HELPER_H

class NormalHelper : public vws::UnivariateHelper<double>
{
public:
	NormalHelper(double z, double lambda2)
	: _z(z), _lambda2(lambda2)
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
	bool s(double x) const {
		return true;
	}
	const NormalHelper& operator=(const NormalHelper& x) {
		_z = x._z;
		_lambda2 = x._lambda2;
		return *this;
	}

private:
	double _z;
	double _lambda2;
};

#endif
