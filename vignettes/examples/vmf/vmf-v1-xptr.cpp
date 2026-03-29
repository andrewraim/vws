// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "texp.h"

typedef vws::FMMProposal<double, vws::RealConstRegion> t_proposal;
typedef Rcpp::XPtr<t_proposal> t_proposal_xptr;
typedef Rcpp::XPtr<vws::UnivariateHelper> t_univariate_helper_xptr;
typedef Rcpp::XPtr<vws::RealConstRegion> t_real_const_region_xptr;

// [[Rcpp::export]]
t_proposal_xptr r_vmf_pre_v1_xptr(double kappa, double d)
{
	/*
    vws::dfdb w = [=](double x, bool log = true)
    {
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
	*/

    // auto w_ptr = std::make_unique<vws::dfdb>(std::move(w));
    auto w_ptr = std::make_unique<vws::dfdb>(
    [=](double x, bool log = true)
    {
        double out = R_NegInf;
        if (std::fabs(x) < 1){
            out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2));
        }
        return log ? out : std::exp(out);
    });
    // auto df_ptr = std::make_unique<fntl::density>(std::move(df));
    // auto pf_ptr = std::make_unique<fntl::cdf>(std::move(pf));
    // auto qf_ptr = std::make_unique<fntl::quantile>(std::move(qf));

    auto df_ptr = std::make_unique<fntl::density>(
    [=](double x, bool log = false) {
        return d_texp(x, kappa, -1, 1, log);
    });

    auto pf_ptr = std::make_unique<fntl::cdf>(
    [=](double q, bool lower = true, bool log = false) {
        return p_texp(q, kappa, -1, 1, lower, log);
    });

    auto qf_ptr = std::make_unique<fntl::quantile>(
    [=](double p, bool lower = true, bool log = false) {
        return q_texp(p, kappa, -1, 1, lower, log);
    });


    // I think we have an issue here because the variables above are defined on
    // the stack ... they won't exist when this function returns.

    // auto* helper_ptr = new vws::UnivariateHelper(df, pf, qf);
    // Rcpp::XPtr<vws::UnivariateHelper> helper = Rcpp::XPtr<vws::UnivariateHelper>(helper_ptr);

	auto helper_ptr = std::make_unique<vws::UnivariateHelper>(*df_ptr, *pf_ptr, *qf_ptr);

    // auto* supp_ptr = new vws::RealConstRegion(-1, 1, w, *helper);
    // auto supp = new Rcpp::XPtr<vws::RealConstRegion>(supp_ptr);

    auto supp_ptr = std::make_unique<vws::RealConstRegion>(-1, 1, *w_ptr, *helper_ptr);

    return t_proposal_xptr(new t_proposal(*supp_ptr), true);
}

// [[Rcpp::export]]
Rcpp::NumericVector refine(t_proposal_xptr h, unsigned int N, double tol = 0)
{
	Rprintf("Made it inside refine\n");
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



/*
#include <Rcpp.h>
#include <vector>

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<int>> create_vector()
{
	std::vector<int>* vec = new std::vector<int>;
	vec->push_back(1);
	vec->push_back(2);
	return Rcpp::XPtr<std::vector<int>>(vec);
}

// [[Rcpp::export]]
void add_value(Rcpp::XPtr<std::vector<int>> ptr, int val)
{
	// Use the pointer to modify the object
	ptr->push_back(val);
}

// [[Rcpp::export]]
void print_vector(Rcpp::XPtr<std::vector<int>> ptr)
{
	for(int i : *ptr) {
		Rcpp::Rcout << i << " ";
	}
	Rcpp::Rcout << std::endl;
}

*/
