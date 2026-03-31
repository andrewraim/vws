// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "texp.h"

/*
* Define a typedef for the proposal
*/
typedef vws::FMMProposal<double, vws::RealConstRegion> t_proposal;

/*
* Define a subclass of FMMProposal that we can expose to R via Modules
*/
class RcppProposal : t_proposal
{
public:
	RcppProposal(double kappa, double d)
	: t_proposal(supp(kappa, d))
	{
	}

	Rcpp::List draw(unsigned int n, unsigned int max_rejects = 10000,
		unsigned int report = 10000);

	Rcpp::NumericVector refine(unsigned int N, double tol = 0) {
		return t_proposal::refine(N, tol);
	}

private:
	vws::RealConstRegion supp(double kappa, double d);
};

/*
* Define module for RcppProposal
*/
RCPP_MODULE(vws_module) {
	Rcpp::class_<RcppProposal>("RcppProposal")
	.constructor<double,double>()
	.method("draw", &RcppProposal::draw)
	.method("refine", &RcppProposal::refine)
	;
}

/*
* Implementation of member functions is below
*/

inline Rcpp::List RcppProposal::draw(
	unsigned int n,
	unsigned int max_rejects,
	unsigned int report)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    auto out = vws::rejection(*this, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects
    );
}

inline vws::RealConstRegion RcppProposal::supp(double kappa, double d)
{
    vws::dfdb w =
    [=](double x, bool log = true) {
        double out = R_NegInf;
        if (std::fabs(x) < 1){
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2));
        }
        return log ? out : std::exp(out);
    };

    fntl::density df = [=](double x, bool log = false) {
        return d_texp(x, kappa, -1, 1, log);
    };

    fntl::cdf pf = [=](double q, bool lower = true, bool log = false) {
        return p_texp(q, kappa, -1, 1, lower, log);
    };

    fntl::quantile qf = [=](double p, bool lower = true, bool log = false) {
        return q_texp(p, kappa, -1, 1, lower, log);
    };

	vws::UnivariateHelper helper(df, pf, qf);
    vws::RealConstRegion out(-1, 1, w, helper);
    return out;
}
