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
		return R::dnorm(x, _z, std::sqrt(_lambda2), log);
	}
	double p(double q, bool lower = true, bool log = false) const {
		return R::pnorm(q, _z, std::sqrt(_lambda2), lower, log);
	}
	double q(double p, bool lower = true, bool log = false) const {
		return R::qnorm(p, _z, std::sqrt(_lambda2), lower, log);
	}
	double s(double x) const {
		return 1.0;
	}
	double w(double x, bool log = false) const {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - _mu, 2.0) / (2*_sigma2);
		}
		return log ? out : exp(out);
	}
	const MyHelper& operator=(const MyHelper& x) {
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
	printf("control.max_rejects = %d\n", control.get_max_rejects());
	printf("control.report_period = %d\n", control.get_report_period());
	printf("control.max_rejects_action = %d\n", control.get_max_rejects_action());

	printf("About to build helper\n");

	MyHelper helper(mu, sigma2, z, lambda2);
	vws::UnivariateConstRegion supp(0.0, R_PosInf, helper);

	printf("About to construct supp\n");

	std::vector<vws::UnivariateConstRegion> regions;
	regions.push_back(supp);

	printf("About to make proposal\n");

	vws::FMMProposal<double, vws::UnivariateConstRegion> h(regions);

	// printf("About to adapt\n");
	// h.adapt(9);

	h.print();

	printf("About to call rejection\n");

	const std::pair<std::vector<double>, std::vector<unsigned int>>& out =
		vws::rejection(h, n, control);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.first,
		Rcpp::Named("rejections") = out.second
	);
}
