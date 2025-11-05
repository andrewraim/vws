# Introduction

`vws` is an R package to support rejection sampling using vertical weighted
strips ([arxiv:2401.09696][arxiv]). Construction of the proposal
distribution and rejection sampling are carried out in C++; sampling functions
may be exposed in R via Rcpp for use in applications. Programming in C++ is 
facilitated using the [fntl][fntl] R package.

[arxiv]: https://arxiv.org/abs/2401.09696
[fntl]: https://cran.r-project.org/package=fntl

See the included vignette for a more in-depth discussion of the package and an
API guide.

# Installation

The `vws` package may be installed directly from Github using a standard R
command like the following.

```r
devtools::install_github("andrewraim/vws", ref = "v0.3.0")
```

Here, `v0.3.0` represents a tagged release; replace it with a later version if
one exists.


# Getting Started

The following example from the vignette generates variates from the density

$$
f(y \mid z, \mu, \sigma^2)
%
\propto
\underbrace{\frac{1}{\lambda \sqrt{2\pi}} \exp\left[
-\frac{1}{2\lambda^2} (z - y)^2
\right]}_{g(y)} \cdot
\underbrace{\frac{1}{y} \exp\left[
-\frac{1}{2\sigma^2} (\log y - \mu)^2
\right] \mathrm{I}(y > 0)}_{w(y)}.
$$

Create the file `example.cpp` with the following contents.

```cpp
// [[Rcpp::depends(vws, fntl)]]
#include "vws.h"

// [[Rcpp::export]]
Rcpp::List example(unsigned int n, double z, double mu,
	double sigma, double lambda, unsigned int N, double tol = 0,
	unsigned int max_rejects = 10000, unsigned int report = 10000)
{
	vws::rejection_args args;
	args.max_rejects = max_rejects;
	args.report = report;
	args.action = fntl::error_action::STOP;

	const vws::dfdb& w =
	[&](double x, bool log = true) {
		double out = R_NegInf;
		if (x > 0) {
			out = -std::log(x) - std::pow(std::log(x) - mu, 2) / (2*sigma*sigma);
		}
		return log ? out : std::exp(out);
	};

	fntl::density df = [&](double x, bool log = false) {
		return R::dnorm(x, z, lambda, log);
	};
	fntl::cdf pf = [&](double q, bool lower = true, bool log = false) {
		return R::pnorm(q, z, lambda, lower, log);
	};
	fntl::quantile qf = [&](double p, bool lower = true, bool log = false) {
		return R::qnorm(p, z, lambda, lower, log);
	};

	vws::UnivariateHelper helper(df, pf, qf);
	vws::RealConstRegion supp(0, R_PosInf, w, helper);
	vws::FMMProposal<double, vws::RealConstRegion> h(supp);

	auto lbdd = h.refine(N - 1, tol);
	const vws::rejection_result<double>& out = vws::rejection(h, n, args);

	return Rcpp::List::create(
		Rcpp::Named("draws") = out.draws,
		Rcpp::Named("rejects") = out.rejects,
		Rcpp::Named("lbdd") = lbdd
	);
}
```

The `example` function may be called through R as follows.

```r
R> Rcpp::sourceCpp("example.cpp")
R> mu = 5; sigma = sqrt(0.5); lambda = 10; y_true = 58; z = 63
R> out = example(n = 1000, z, mu, sigma, lambda, N = 50, tol = 0.10)
R> head(out$draws)
```
