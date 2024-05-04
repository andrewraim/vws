// [[Rcpp::depends(vws)]]
#include "vws.h"

class MyHelper : public vws::UnivariateHelper<double>
{
public:
	MyHelper(double mu, double sigma2, double lambda2)
		: _mu(mu), _sigma2(sigma2), _lambda2(lambda2)
	{
	}

	double d(double x, bool log = false) const {
		return 0;
	}
	double p(double q, bool lower = true, bool log = false) const {
		return 0;
	}
	double q(double p, bool lower = true, bool log = false) const {
		return 0;
	}
	double s(double x) const {
		return 0;
	}
	double w(double x, bool log = false) const {
		double out = x > 0 ? -std::log(x) - std::pow(std::log(x) - _mu, 2.0) / (2*_sigma2) + std::log(x > 0) : R_NegInf;
		if (log) { return out; } else { return exp(out); }
	}

private:
	double _mu;
	double _sigma2;
	double _lambda2;
};

// [[Rcpp::export]]
int hello(double x)
{
	printf("Hello world %g\n", x);

	vws::RejectionControl control;
	printf("control.max_rejects_action = %d\n", control.get_max_rejects_action());

	MyHelper helper(0.0, 1.0, 5.0);
	vws::UnivariateConstRegion supp(0, R_PosInf, helper);

	return 0;
}
