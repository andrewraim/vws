// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "LinearVWSRegion.h"

// [[Rcpp::export]]
Rcpp::List r_bessel_v3(unsigned int n, double lambda, double nu, unsigned int N,
	double lo, double hi, double tol = 0, unsigned int max_rejects = 10000,
	unsigned int report = 10000,
	const Rcpp::NumericVector& x = Rcpp::NumericVector::create())
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = report;

	LinearVWSRegion supp(lo, hi, lambda, nu);
	vws::FMMProposal<double, LinearVWSRegion> h(supp);

	auto lbdd = h.refine(N - 1, tol);
	auto out = vws::rejection(h, n, args);

	h.print(N);

    // Evaluate the majorized weight function, minorized weight function, and
    // proposal on the given values of x. The caller can plot these.
	unsigned int m = x.size();
    Rcpp::NumericVector hx(m);
    Rcpp::NumericVector wmajx(m);
    for (unsigned int i = 0; i < m; i++) {
        hx(i) = h.d(x(i));
        wmajx(i) = h.w_major(x(i));
    }


    /*
    * Get the bounds and coefficients for the linear majorizer. The caller can
    * use these to plot the majorizer.
    */
    Rcpp::NumericVector lower(N);
    Rcpp::NumericVector upper(N);
    Rcpp::NumericVector beta0_max(N);
    Rcpp::NumericVector beta1_max(N);
    Rcpp::NumericVector beta0_min(N);
    Rcpp::NumericVector beta1_min(N);

    std::set<LinearVWSRegion>::const_iterator itr = h.regions_begin();
    unsigned int i = 0;
    for (; itr != h.regions_end(); itr++) {
		lower(i) = itr->get_lower();
    	upper(i) = itr->get_upper();
    	beta0_max(i) = itr->get_beta0_max();
    	beta1_max(i) = itr->get_beta1_max();
    	beta0_min(i) = itr->get_beta0_min();
    	beta1_min(i) = itr->get_beta1_min();
    	i++;
    }

    const Rcpp::DataFrame& df_weight = Rcpp::DataFrame::create(
    	Rcpp::Named("lo") = lower,
    	Rcpp::Named("hi") = upper,
    	Rcpp::Named("beta0_max") = beta0_max,
    	Rcpp::Named("beta1_max") = beta1_max,
    	Rcpp::Named("beta0_min") = beta0_min,
    	Rcpp::Named("beta1_min") = beta1_min
    );

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = lbdd,
        Rcpp::Named("hx") = hx,
        Rcpp::Named("wmajx") = wmajx,
        Rcpp::Named("df_weight") = df_weight
	);
}
