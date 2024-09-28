#ifndef VWS_RECT_H
#define VWS_RECT_H

#include <Rcpp.h>
#include "logit.h"

namespace vws {

inline double rect(double z, double a, double b)
{
	double x;

	if (std::isinf(a) && std::isinf(b) && a < 0 && b > 0) {
		x = z;
	} else if (std::isinf(a) && a < 0) {
		x = logit(z / b);
	} else if (std::isinf(b) && b > 0) {
		x = log(z - a);
	} else if (std::isfinite(a) && std::isfinite(b)) {
		x = logit(z / (b - a) - a);
	} else {
		Rcpp::stop("Invalid endpoints");
	}

	return x;
}


inline double inv_rect(double x, double a, double b)
{
	double z;

	if (std::isinf(a) && std::isinf(b) && a < 0 && b > 0) {
		z = x;
	} else if (std::isinf(a) && a < 0) {
		z = b * inv_logit(x);
	} else if (std::isinf(b) && b > 0) {
		z = exp(x) + a;
	} else if (std::isfinite(a) && std::isfinite(b)) {
		z = (b - a) * inv_logit(x) + a;
	} else {
		Rcpp::stop("Invalid endpoints");
	}

	return z;
}

inline Rcpp::NumericVector rect(const Rcpp::NumericVector& z,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	unsigned int n = z.size();
	if (n != a.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (n != b.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (Rcpp::is_true(Rcpp::any(a > b))) {
		Rcpp::stop("any(a > b)");
	}

	Rcpp::NumericVector x(n);

	for (unsigned int i = 0; i < n; i++) {
		if (std::isinf(a(i)) && std::isinf(b(i)) && a(i) < 0 && b(i) > 0) {
			x(i) = z(i);
		} else if (std::isinf(a(i)) && a(i) < 0) {
			x(i) = logit(z(i) / b(i));
		} else if (std::isinf(b(i)) && b(i) > 0) {
			x(i) = log(z(i) - a(i));
		} else if (std::isfinite(a(i)) && std::isfinite(b(i))) {
			x[i] = logit(z(i) / (b(i) - a(i)) - a(i));
		} else {
			Rcpp::stop("Invalid endpoints");
		}
	}

	return x;
}

inline Rcpp::NumericVector inv_rect(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& a, const Rcpp::NumericVector& b)
{
	unsigned int n = x.size();
	if (n != a.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (n != b.size()) {
		Rcpp::stop("Dimension mismatch");
	}
	if (is_true(Rcpp::any(a > b))) {
		Rcpp::stop("any(a > b)");
	}

	Rcpp::NumericVector z(n);

	for (unsigned int i = 0; i < n; i++) {
		if (std::isinf(a(i)) && std::isinf(b(i)) && a(i) < 0 && b(i) > 0) {
			z(i) = x(i);
		} else if (std::isinf(a(i)) && a(i) < 0) {
			z(i) = b(i) * inv_logit(x(i));
		} else if (std::isinf(b(i)) && b(i) > 0) {
			z(i) = exp(x(i)) + a(i);
		} else if (std::isfinite(a(i)) && std::isfinite(b(i))) {
			z(i) = (b(i) - a(i)) * inv_logit(x(i)) + a(i);
		} else {
			Rcpp::stop("Invalid endpoints");
		}
	}

	return z;
}

}

#endif
