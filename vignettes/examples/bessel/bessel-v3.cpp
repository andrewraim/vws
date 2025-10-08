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

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = lbdd,
        Rcpp::Named("hx") = hx,
        Rcpp::Named("wmajx") = wmajx
	);
}
