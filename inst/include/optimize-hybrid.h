#ifndef VWS_OPTIMIZE_HYBRID_H
#define VWS_OPTIMIZE_HYBRID_H

#include <Rcpp.h>
#include "fntl.h"
#include "logit.h"

namespace vws {

struct optimize_hybrid_result {
	double par;
	double value;
	std::string method;
};

/*
* A hybrid optimization for univariate problems. Use Brent if both bounds are
* finite; otherwise, use L-BFGS-B. In the latter cause, if one of the bounds is
* finite, enforce it with a transformation.
*/
inline optimize_hybrid_result optimize_hybrid(const fntl::dfd& f, double init,
	double lower, double upper, bool maximize, unsigned maxiter = 100000)
{
	if (lower > upper) {
		Rcpp::stop("lower > upper");
	}

	double f_lower = f(lower);
	double f_upper = f(upper);
	bool f_lower_pos_inf = (f_lower > 0) && std::isinf(f_lower);
	bool f_lower_neg_inf = (f_lower < 0) && std::isinf(f_lower);
	bool f_upper_pos_inf = (f_upper > 0) && std::isinf(f_upper);
	bool f_upper_neg_inf = (f_upper < 0) && std::isinf(f_upper);

	optimize_hybrid_result out;

	if (maximize && f_lower_pos_inf) {
		out.par = lower;
		out.value = R_PosInf;
		out.method = "Lower Limit Inf";
		return out;
	}

	if (maximize && f_upper_pos_inf) {
		out.par = upper;
		out.value = R_PosInf;
		out.method = "Upper Limit Inf";
		return out;
	}

	if (!maximize && f_lower_neg_inf) {
		out.par = lower;
		out.value = R_NegInf;
		out.method = "Lower Limit NegInf";
		return out;
	}

	if (!maximize && f_upper_neg_inf) {
		out.par = upper;
		out.value = R_NegInf;
		out.method = "Upper Limit NegInf";
		return out;
	}

	if (std::isfinite(lower) && std::isfinite(upper)) {
		// If both endpoints are finite, use Brent's method
	    fntl::optimize_args args;
	    args.fnscale = maximize ? -1 : 1;
	    args.maxiter = maxiter;

		const auto& opt_out = fntl::optimize_brent(f, lower, upper, args);
		out.par = opt_out.par;
		out.value = opt_out.value;
		out.method = "Brent";
	} else {
		// If one or both endpoints are infinite, use L-BFGS-B method
		fntl::lbfgsb_args args;
		args.fnscale = maximize ? -1 : 1;
		args.maxit = maxiter;

		// Transform to the interval (a,b]
		const fntl::dfd& tx = [&](double x) {
			if (std::isinf(lower) && std::isinf(upper) && lower < 0 && upper > 0) {
				return x;
			} else if (std::isinf(lower) && lower < 0) {
				return upper * logit(x);
			} else if (std::isinf(upper) && upper > 0) {
				return std::exp(x) + lower;
			} else {
				return (upper - lower) * logit(x) + lower;
			}
		};

		// Compose function with transformation
	    const fntl::dfv& ff = [&](const Rcpp::NumericVector& x) {
			return f(tx(x(0)));
		};

		const auto& vinit = Rcpp::NumericVector::create(init);
		const fntl::lbfgsb_result& opt_out = fntl::lbfgsb(vinit, ff, args);

		if (opt_out.status != fntl::lbfgsb_status::OK) {
			Rcpp::warning("L-BFGS-B: convergence status was %d",
				fntl::to_underlying(opt_out.status));
		}

		out.par = tx(opt_out.par[0]);
		out.value = opt_out.value;
		out.method = "L-BFGS-B";
	}

	// In case the function is strictly increasing or decreasing, check the
	// objective value at the endpoints.
	if (maximize) {
		if (f_lower > out.value) {
			out.par = lower;
			out.value = f_lower;
			out.method = "Max at Lower Limit";
		}
		if (f_upper > out.value) {
			out.par = upper;
			out.value = f_upper;
			out.method = "Max at Upper Limit";
		}
	} else {
		if (f_lower < out.value) {
			out.par = lower;
			out.value = f_lower;
			out.method = "Min at Lower Limit";
		}
		if (f_upper < out.value) {
			out.par = upper;
			out.value = f_upper;
			out.method = "Min at Upper Limit";
		}
	}

	return out;
}

}

#endif
