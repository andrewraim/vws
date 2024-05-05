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
		return R::dnorm(x, _z, sqrt(_lambda2), log);
	}
	double p(double q, bool lower = true, bool log = false) const {
		return R::pnorm(q, _z, sqrt(_lambda2), lower, log);
	}
	double q(double p, bool lower = true, bool log = false) const {
		return R::qnorm(p, _z, sqrt(_lambda2), lower, log);
	}
	double s(double x) const {
		return 1.0;
	}
	double w(double x, bool log = false) const {
		double out = x > 0 ? -std::log(x) - std::pow(std::log(x) - _mu, 2.0) / (2*_sigma2) + std::log(x > 0) : R_NegInf;
		if (log) { return out; } else { return exp(out); }
	}
	UnivariateHelper<T> operator=(const UnivariateHelper<T>& x) {
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

// [[Rcpp::export]]
Rcpp::List r_lognormal_normal(unsigned int n, double z, double mu, double sigma2,
	double lambda2)
{
	vws::RejectionControl control;
	printf("control.max_rejects_action = %d\n", control.get_max_rejects_action());

	MyHelper helper(0.0, 1.0, -10, 5.0);
	vws::Region<double> supp(0, R_PosInf, helper);

	std::vector<vws::Region<double>> regions;
	regions.push_back(supp);

	vws::FMMProposal<double> h(regions);

	const std::pair<std::vector<double>, Rcpp::IntegerVector> out = vws::rejection(h, n, control);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejections") = out.second
	);
}
