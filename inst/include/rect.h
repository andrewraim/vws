#ifndef VWS_RECT_H
#define VWS_RECT_H

#include <Rcpp.h>
#include "logit.h"

namespace vws {

Rcpp::NumericVector rect(const Rcpp::NumericVector& z,
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

Rcpp::NumericVector inv_rect(const Rcpp::NumericVector& x,
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
