// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "texp.h"

/*
* Utility functions to set and check tags for an XPtr
*/

template <typename T>
void assert_tag(Rcpp::XPtr<T> x, const char* tag)
{
	if (Rf_isNull(R_ExternalPtrTag(x))) {
		Rcpp::stop("XPtr does not have a tag");
	}

	if (Rcpp::as<Rcpp::String>(R_ExternalPtrTag(x)) != Rcpp::String(tag)) {
		Rcpp::stop("XPtr is not tagged as %s", tag);
	}
}

template <typename T>
void set_tag(Rcpp::XPtr<T> x, const char* tag)
{
	R_SetExternalPtrTag(x, Rcpp::wrap(Rcpp::String(tag)));
}

/*
* Define typedefs for the proposal and its XPtr
*/
typedef vws::FMMProposal<double, vws::RealConstRegion> t_proposal;
typedef Rcpp::XPtr<t_proposal> t_proposal_xptr;

// [[Rcpp::export]]
t_proposal_xptr vmf_pre_v1_xptr(double kappa, double d)
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

    auto p = new t_proposal(supp);
    auto out = t_proposal_xptr(p, true);
    set_tag(out, "t_proposal");
	return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector refine(t_proposal_xptr h, unsigned int N, double tol = 0)
{
	assert_tag(h, "t_proposal");
    return h->refine(N, tol);
}

// [[Rcpp::export]]
Rcpp::List draw(t_proposal_xptr h, unsigned int n,
	unsigned int max_rejects = 10000, unsigned int report = 10000)
{
	assert_tag(h, "t_proposal");

    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    auto out = vws::rejection(*h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects
    );
}

/*
* This function returns another XPtr which is incompatible with t_proposal, for
* testing.
*/

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double>> test_vector()
{
	std::vector<double>* p = new std::vector<double>();
	p->push_back(10);
	p->push_back(10);
	p->push_back(12);
	return Rcpp::XPtr<std::vector<double>>(p);
}

