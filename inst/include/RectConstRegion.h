#ifndef VWS_RECT_CONST_REGION_H
#define VWS_RECT_CONST_REGION_H

#include <Rcpp.h>
#include <memory>
#include "fntl.h"
#include "Region.h"
#include "UnivariateHelper.h"

namespace vws {

//' Rectangular Region with Constant Majorizer
//'
//' A class which represents a region based on rectagular intervals with a
//' constant majorizer for the weight function. This version is for continuous
//' supports.
//'
//' @field w Weight function for the target distribution. Its expected interface
//' is \code{w(x, log = TRUE)} so that results are returned on the log-scale by
//' default.
//'
//' @examples
//' # Define base distribution and weight function
//' g = normal_univariate_helper(mean = 0, sd = 5)
//' w = function(x, log = FALSE) { dlnorm(10 - x, meanlog = 5, sdlog = 2, log) }
//'
//' reg = RectConstRegion$new(-Inf, 10, w, g)
//' print(reg)
//'
//' out = reg$bifurcate(0)
//' print(out[[1]])
//' print(out[[2]])
//'
//' @export
class UnivariateConstRegion : public Region<double>
{
protected:
	Rcpp::NumericVector _a;
	Rcpp::NumericVector _b;
	const mv_weight_function* _w;
	const std::vector<UnivariateHelper<double>>* _helpers;
	double _log_w_max;
	double _log_w_min;
	double _log_prob;

public:
	//' @param a Lower and upper limit of interval.
	//' @param w Weight function for the target distribution.
	//' @param g An object created by \code{univariate_helper}.
	RectConstRegion(const Rcpp::NumericVector& a,
		const mv_weight_function& w,
		const std::vector<UnivariateHelper<double>>& helpers);

	//' @param a Lower limit of interval.
	//' @param b Upper limit of interval.
	//' @param w Weight function for the target distribution.
	//' @param g An object created by \code{univariate_helper}.
	UnivariateConstRegion(const Rcpp::NumericVector& a,
		const Rcpp::NumericVector& b, const mv_weight_function& w,
		const std::vector<UnivariateHelper<double>>& helpers);

	//' @description
	//' Density function \eqn{g} for the base distribution.
	//' @param x Density argument.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double d_base(const Rcpp::NumericVector& x, bool log = false) const;

	//' @description
	//' Generate a draw from \eqn{g_j} specific to this region.
	//' @param n Number of draws to generate.
	//' @return A list of draws, with one draw per list element.
	std::vector<Rcpp::NumericVector> r(unsigned int n) const;

	//' @description
	//' Density of \eqn{g_j} specific to this region.
	//' @param x Density argument.
	double d(const Rcpp::NumericVector& x, bool log = false) const;

	//' @description
	//' Test if given \code{x} is in the support for the \eqn{g_j} specific to
	//' this region.
	//' @param x Density argument.
	bool s(const Rcpp::NumericVector& x) const;

	double w(const Rcpp::NumericVector& x, bool log = true) const;

	//' @description
	//' Return a logical value indicating whether this region is bifurcatable.
	bool is_bifurcatable() const;

	//' @description
	//' The quantity \eqn{\overline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double get_xi_upper(bool log = true) const;

	//' @description
	//' The quantity \eqn{\underline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double get_xi_lower(bool log = true) const;

	//' @description
	//' A string that describes the region.
	std::string description() const;

	//' @description
	//' Print a description of the region.
	void print() const;

	//' @description
	//' Maximize or minimize the function \eqn{w(x)} over this region. Optimization
	//' is carried out with \code{optim} using arguments
	//' \code{method = "Nelder-Mead"} and
	//' \code{control = list(maxit = 100000, warn.1d.NelderMead = FALSE)}.
	//' @param maximize logical; if \code{TRUE} do maximization. Otherwise do
	//' minimization.
	//' @param log logical; if \code{TRUE} return optimized value of \eqn{\log w(x)}.
	//' Otherwise return optimized value of \eqn{w(x)}.
	double optimize(bool maximize = true, bool log = true) const;

	//' @description
	//' Majorized weight function \eqn{\overline{w}_j} for this region.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double w_major(const Rcpp::NumericVector& x, bool log = true) const;

	Rcpp::NumericVector midpoint() const;

	std::pair<RectConstRegion,RectConstRegion> bifurcate() const;

	//' @description
	//' Bifurcate this region into two regions. Use \code{x} as the bifurcation
	//' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
	//' @param x An optional bifurcation point.
	std::pair<RectConstRegion,RectConstRegion> bifurcate(unsigned int l, double x) const;

	RectConstRegion singleton(const Rcpp::NumericVector& x) const;

	bool operator<(const RectConstRegion& x) const
	{
		for (unsigned int l = 0; l < k; l++) {
			if (_a(l) < x._a(l) && _b(l) > x._a(l)) {
				Rcpp::stop("Overlapping regions found");
			}
			if (x._a(l) < _a(l) && x._b(l) > _a(l)) {
				Rcpp::stop("Overlapping regions found");
			}

			if (_a(l) < x._a(l)) {
				return true;
			} else if (_a(l) > x._a(l) {
				return false;
			}
		}

		return false;
	}

	bool operator==(const RectConstRegion& x) const
	{
		return Rcpp::all(_a == x._a) && Rcpp::all(_b == x._b);
	}

	const RectConstRegion& operator=(const RectConstRegion& x)
	{
		_a = x._a;
		_b = x._b;
		_w = x._w;
		_helpers = x._helpers;
		_log_w_max = x._log_w_max;
		_log_w_min = x._log_w_min;
		_log_prob = x._log_prob;
		return *this;
	}
};

RectConstRegion::RectConstRegion(const Rcpp::NumericVector& a,
	const mv_weight_function& w, const std::vector<UnivariateHelper<double>>& helpers)
: _a(a), _b(a), _w(&w), _helpers(&helpers)
{
	_log_w_max = (*_w)(a, true);
	_log_w_min = (*_w)(a, true);
	_log_prob = R_NegInf;
}

RectConstRegion::RectConstRegion(const Rcpp::NumericVector& a,
	const Rcpp::NumericVector& b, const mv_weight_function& w,
	const std::vector<UnivariateHelper<double>>& helpers)
: _a(a), _b(b), _w(&w), _helpers(&helpers)
{
	unsigned int k = helpers.size();

	if (k != length(a) || k == length(b)) {
		Rcpp::stop("Dimension mismatch");
	}

	if (Rcpp::any(a > b)) {
		Rcpp::stop("a > b");
	}

	_log_w_max = optimize(true);
	_log_w_min = optimize(false);

	// Compute Prob((a,b]) under dist g on the log scale
	_log_prob = 0;
	for (unsigned int l = 0; l < k; i++) {
		lpa = g[l]->p(a(l), true, true);
		lpb = g[l]->p(b(l), true, true);
		_log_prob += log_sub2_exp(lpb, lpa);
	}
}

double RectConstRegion::d_base(const Rcpp::NumericVector& x, bool log) const
{
	unsigned int k = length(_g);
	if (length(x) != k) {
		Rcpp::stop("Dimension mismatch");
	}

	out = 0;
	for (unsigned int l = 0; l < k; l++) {
		out += _helpers[l]->d(x(l), true);
	}

	return log ? out : exp(out);
}

std::vector<Rcpp::NumericVector> RectConstRegion::r(unsigned int n) const
{
	unsigned int k = length(_helpers);
	Rcpp::NumericMatrix x(n, k);

	for (unsigned int l = 0; l < k; l++) {
		const Rcpp::NumericVector& u = Rcpp::runif(n);
		double lpa = _helper[l]->p(_a(l), true);
		double lpb = _helper[l]->p(_b(l), true);
		double log_prob_l = log_sub2_exp(lpb, lpa);
		const Rcpp::NumericVector& log_p = log_add2_exp(log_prob_l + Rcpp::log(u), rep(lpa, n));
		x.column(l) = _helper[l]->q(log_p, true);
	}

	std::vector<Rcpp::NumericVector> out;
	for (unsigned int i = 0; i < k; i++) {
		out.push_back(x.row(i));
	}

	return out;
}

double RectConstRegion::d(const Rcpp::NumericVector& x, bool log) const
{
	unsigned int k = _helpers.size();

	if (k != length(x)) {
		Rcpp::stop("Dimension mismatch");
	}

	double out = 0;

	for (unsigned int l = 0; l < k; l++) {
		if (a(l) < x(l) & x(l) <= b(l) & _helpers[l]->s(x(l))) {
			double lpa = _helpers[l]->p(a(l), true);
			double lpb = _helpers[l]->p(b(l), true);
			double lpx = _helpers[l]->d(x(l), true);
			out += lpx - log_sub2_exp(lpb, lpa);
		}
	}

	return log ? out : exp(out);
}

bool RectConstRegion::s(const Rcpp::NumericVector& x) const
{
	unsigned int k = _helpers.size();

	if (k != length(x)) {
		Rcpp::stop("Dimension mismatch");
	}

	bool out = true;

	for (unsigned int l = 0; l < k; l++) {
		out &= (_a(l) < x(l) && x(l) <= _b(l)) & _helper[l]->s(x(l));
	}

	return out;
}

double RectConstRegion::w(const Rcpp::NumericVector& x, bool log) const
{
	return (*_w)(x, log);
}

double RectConstRegion::w_major(const Rcpp::NumericVector& x, bool log) const
{
	double out = _helper->s(x) ? _log_w_max : R_NegInf;
	return log ? out : exp(out);
}

Rcpp::NumericVector RectConstRegion::midpoint() const
{
	Rcpp::NumericVector out(k);

	for (unsigned int l = 0; l < k; l++) {
		if (std::isinf(_a(k)) && std::isinf(_b(k)) && _a(k) < 0 && _b(k) > 0) {
			// In this case, we have an interval (-Inf, Inf). Make a split at zero.
			out(k) = 0;
		} else if (std::isinf(_a(k)) && _a(k) < 0) {
			// Left endpoint is -Inf. Split based on right endpoint.
			out(k) = _b(k) - abs(_b) - 1;
		} else if (std::isinf(_b(k)) && _b(k) > 0) {
			// Right endpoint is Inf. Split based on left endpoint.
			out(k) = _a(k) + std::fabs(_a(k)) + 1;
		} else {
			out(k) = (_a(k) + _b(k)) / 2;
		}
	}

	return out;
}

std::pair<RectConstRegion,RectConstRegion>
RectConstRegion::bifurcate() const
{
	unsigned int k = _helpers.size();
	const Rcpp::NumericVector& cuts = midpoint();

	if (true) {
		std::vector<RectConstRegion> pairs_left;
		std::vector<RectConstRegion> pairs_right;
		Rcpp::NumericVector reduction(k);

		for (unsigned int l = 0; l < k, l++) {
		 	Rcpp::NumericVector a1 = _a;
		 	Rcpp::NumericVector a2 = _a;
		 	Rcpp::NumericVector b1 = _b;
		 	Rcpp::NumericVector b2 = _b;
		 	a2(l) = cuts(l);
		 	b1(l) = cuts(l);

			RectConstRegion r1(a1, b1, _w, _helpers);
			RectConstRegion r2(a2, b2, _w, _helpers);
			pairs_left.push_back(r1);
			pairs_right.push_back(r2);

			double lp0 = log_sub2_exp(self$xi_upper(), xi_lower());
			double lp1 = log_sub2_exp(r1.xi_upper(), r1.xi_lower());
			double lp2 = log_sub2_exp(r2.xi_upper(), r2.xi_lower());

			if (lp0 > log_add2_exp(lp1, lp2)) {
				reduction(l) = log_sub2_exp(lp0, log_add2_exp(lp1, lp2));
			} else {
				reduction(l) = R_NegInf;
			}
		}

		unsigned int idx = r_categ(1, reduction, true);
		return std::make_pair(pairs_left[idx], pairs_right[idx]);
	} else {
		// Sample a cut orientation randomly
		Rcpp::NumericVector p(k);
		p.fill(1.0 / k);
		unsigned int idx = r_categ(1, p);
		Rcpp::NumericVector a1 = _a;
		Rcpp::NumericVector a2 = _a;
		Rcpp::NumericVector b1 = _b;
		Rcpp::NumericVector b2 = _b;
		a2(idx) = cuts(idx);
		b1(idx) = cuts(idx);
		RectConstRegion r1(a1, b1, _w, _helpers);
		RectConstRegion r1(a2, b2, _w, _helpers);
		return std::make_pair(r1, r2);
	}
}

std::pair<RectConstRegion,RectConstRegion>
RectConstRegion::bifurcate(unsigned int l, double x) const
{
	a1 = _a;
	a2 = _a;
	b1 = _b;
	b2 = _b;
	a2[l] = x;
	b1[l] = x;

	RectConstRegion r1(a1, b1, *_w, *_helpers);
	RectConstRegion r2(a2, b2, *_w, *_helpers);
	return std::make_pair(r1, r2);
}

RectConstRegion RectConstRegion::singleton(const Rcpp::NumericVector& x) const
{
	return RectConstRegion(x, *_w, *_helper);
}

bool RectConstRegion::is_bifurcatable() const
{
	return true;
}

double RectConstRegion::get_xi_upper(bool log) const
{
	double log_xi_upper = _log_w_max + _log_prob;
	return log ? log_xi_upper : exp(log_xi_upper);
}

double RectConstRegion::get_xi_lower(bool log) const
{
	double log_xi_lower = _log_w_min + _log_prob;
	return log ? log_xi_lower : exp(log_xi_lower);
}

std::string RectConstRegion::description() const
{
	std::string out;

	for (unsigned int l = 0; l < k; l++) {
		if (l > 0) { out += " x "; }
		out += "(" + _a(l) + "," + _b(l) + "]";
	}

	return out;
}

void RectConstRegion::print() const
{
	std::string out = description();
	printf("Rect Const Region %s\n", out_c_str());
}

double RectConstRegion::optimize(bool maximize, bool log) const
{
	fntl::lbfgsb_args args;
	args.maxit = 100000;
	args.fnscale = maximize ? -1 : 1;

	std::function<double(Rcpp::NumericVector)> f_opt =
	[&](const Rcpp::NumericVector& x) -> double {
		const Rcpp::NumericVector& tx = inv_rect(x, _a, _b);
		return w(tx, true);
	}

	Rcpp::NumericVector init(k);

	auto opt_out& = fntl::lbfgsb(init, f_opt, args);
	if (opt_out.convergence != 0 && opt_out.convergence != 52) {
		Rcpp::warning("opt_out: convergence status was ", opt_out.convergence);
	}
	double out = opt_out.value;
	return log ? out : exp(out);
}

}

#endif
