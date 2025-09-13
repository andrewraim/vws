// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"
#include "texp.h"

// [[Rcpp::export]]
Rcpp::List r_vmf_pre_v2(unsigned int n, double kappa, double d, unsigned int N,
    double tol = 0, unsigned int max_rejects = 10000, unsigned int report = 1000)
{
    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = report;

    const vws::weight_dfd& w =
    [&](double x, bool log = true) {
    	double out = R_NegInf;
    	if (std::fabs(x) < 1){
        	out = 0.5 * (d - 3) * std::log1p(-std::pow(x, 2));
    	}
        return log ? out : std::exp(out);
    };

    fntl::density df = [&](double x, bool log = false) {
        return d_texp(x, kappa, -1, 1, log);
    };
    fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
        return p_texp(q, kappa, -1, 1, lower, log);
    };
    fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
        return q_texp(p, kappa, -1, 1, lower, log);
    };

    const vws::optimizer& maxopt = [&](const vws::weight_dfd& w,
        double lo, double hi, bool log)
    {
    	double out;
    	if (lo <= 0 && 0 < hi) {
    		out = w(0, true);
    	} else if (hi < 0) {
    		out = w(hi, true);
    	} else {
    		out = w(lo, true);
    	}

    	/*
        Rprintf("Begin maxopt\n");
        double A = 1, B = (d-3) / kappa, C = -1;
        double x_root1 = ( -B + std::sqrt(B*B - 4*A*C) ) / (2*A);
        double x_root2 = ( -B - std::sqrt(B*B - 4*A*C) ) / (2*A);
        double x_root = -1 < x_root1 && x_root1 <= 1 ? x_root1 : x_root2;
        if (-1 >= x_root || x_root > 1) { Rcpp::stop("Invalid root"); }

        double out;
        if (x_root < lo) {
            out = w(lo, true);
        } else if (x_root > hi) {
            out = w(hi, true);
        } else {
            out = w(x_root, true);
        }

        Rprintf("End maxopt, out = %g\n", out);
    	*/
        return log ? out : std::exp(out);
    };

    const vws::optimizer& minopt = [&](const vws::weight_dfd& w,
        double lo, double hi, bool log)
    {
    	// Rprintf("Begin minopt\n");
    	double w_lo = w(lo, true);
    	double w_hi = w(hi, true);
    	double out = std::min(w_lo, w_hi);
    	// Rprintf("End minopt, out = %g\n", out);
        return log ? out : std::exp(out);
    };

    vws::UnivariateHelper helper(df, pf, qf);
    // Rprintf("Begin constructing supp\n");
    vws::RealConstRegion supp(-1, 1, w, helper, maxopt, minopt);
    // Rprintf("Begin constructing h\n");
    vws::FMMProposal<double, vws::RealConstRegion> h({ supp });
	// Rprintf("Finish constructing h\n");

    // Rprintf("Begin refining h\n");
    auto lbdd = h.refine(N - 1);
    // Rprintf("Begin rejection\n");
    auto out = vws::rejection(h, n, args);

    return Rcpp::List::create(
        Rcpp::Named("draws") = out.draws,
        Rcpp::Named("rejects") = out.rejects,
        Rcpp::Named("lbdd") = lbdd
    );
}
