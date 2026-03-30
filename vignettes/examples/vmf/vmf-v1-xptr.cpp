// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "texp.h"

// Define types for the proposal and its XPtr
typedef vws::FMMProposal<double, vws::RealConstRegion> t_proposal;
typedef Rcpp::XPtr<t_proposal> t_proposal_xptr;

// [[Rcpp::export]]
t_proposal_xptr r_vmf_pre_v1_xptr(double kappa, double d)
{
    vws::dfdb w =
    [=](double x, bool log = true) {
        double out = R_NegInf;
        if (std::fabs(x) < 1){
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2));
        }
        return log ? out : std::exp(out);
    };

    /*
    * We need to capture by value here (using [=]) because kappa and d are
    * defined locally in this function.
    */

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
    vws::RealConstRegion supp(-1, 1, w, helper);

    auto out = new t_proposal(supp);
    return t_proposal_xptr(out, true);
}

// [[Rcpp::export]]
Rcpp::NumericVector refine(t_proposal_xptr h, unsigned int N, double tol = 0)
{
    return h->refine(N - 1, tol);
}

// [[Rcpp::export]]
Rcpp::List draw(t_proposal_xptr h, unsigned int n,
	unsigned int max_rejects = 10000, unsigned int report = 10000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    auto out = vws::rejection(*h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects
    );
}

// [[Rcpp::export]]
void destruct(t_proposal_xptr h)
{
	/*
	* TBD: We don't have a great way to access the original pointers to delete
	* them right now. We can add accessors to Region or FMMProposal, but it's
	* not good to be able to delete an object's member data.
	*
	* This might be a good time to change the pointers to shared_ptr ...
	*/
	// delete h->regions_begin()->get_w;

	/*
	* TBD: Same issue here. We might be able to have Helper take care of
	* deleting the additional function pointers it has as members ...
	*/
	// delete h->regions_begin()->get_helper;
}
