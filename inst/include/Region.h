#ifndef REGION_H
#define REGION_H

#include <Rcpp.h>

namespace vws {

// A Region contains all of the problem-specific logic for the vws sampler.
// This is an "abstract" R6 class that defines the interface for a Region.
template <class T>
class Region
{
public:
	//' Density function \eqn{g} for the base distribution.
	//' @param x Density argument.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	virtual double d_base(const T& x, bool log = false) const = 0;

	//' @description
	//' Weight function \eqn{w}.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	virtual double w(const T& x, bool log = true) const = 0;

	//' @description
	//' Generate a draw from \eqn{g_j} specific to this region.
	//' @param n Number of draws to generate.
	//' @return A list of draws, with one draw per list element.
	virtual std::vector<T> r(unsigned int n) const = 0;

	//' @description
	//' Density of \eqn{g_j} specific to this region.
	//' @param x Density argument.
	virtual double d(const T& x) const = 0;

	//' @description
	//' Test if given \code{x} is in the support for the \eqn{g_j} specific to this
	//' region.
	//' @param x Density argument.
	virtual bool s(const T& x) const = 0;

	//' @description
	//' Majorized weight function \eqn{\overline{w}_j} for this region.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	virtual double  w_major(const T& x, bool log = true) const = 0;

	//' @description
	//' Bifurcate this region into two regions. Use \code{x} as the bifurcation
	//' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
	//' @param x An optional bifurcation point.
	virtual std::pair<T,T> bifurcate() const = 0;

	virtual std::pair<T,T> bifurcate(const T& x) const = 0;

	//' @description
	//' Return a logical value indicating whether this region is bifurcatable.
	virtual bool is_bifurcatable() const = 0;

	//' @description
	//' The quantity \eqn{\overline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	virtual double xi_upper(bool log = true) const = 0;

	//' @description
	//' The quantity \eqn{\underline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	virtual double xi_lower(bool log = true) const = 0;

	//' @description
	//' A string that describes the region.
	virtual std::string description() const = 0;

	//' @description
	//' Print a description of the region.
	virtual void print() const;
};

}

#endif
