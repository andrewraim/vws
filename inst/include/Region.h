#ifndef REGION_H
#define REGION_H

#include <Rcpp.h>
#include <memory>
#include <type_traits>

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
	double d_base(const T& x, bool log = false) const {
		static_assert(false, "Specialized implementation for template is needed");
		return 0;
	}

	//' @description
	//' Generate a draw from \eqn{g_j} specific to this region.
	//' @param n Number of draws to generate.
	//' @return A list of draws, with one draw per list element.
	std::vector<T> r(unsigned int n) const {
		static_assert(false, "Specialized implementation for template is needed");
		std::vector<T> out;
		return out;
	}

	//' @description
	//' Density of \eqn{g_j} specific to this region.
	//' @param x Density argument.
	double d(const T& x, bool log = false) const {
		static_assert(false, "Specialized implementation for template is needed");
		return 0;
	}

	//' @description
	//' Test if given \code{x} is in the support for the \eqn{g_j} specific to this
	//' region.
	//' @param x Density argument.
	bool s(const T& x) const {
		static_assert(false, "Specialized implementation for template is needed");
		return false;
	}

	double w(const T& x, bool log = true) const {
		static_assert(false, "Specialized implementation for template is needed");
		return 0;
	}

	//' @description
	//' Majorized weight function \eqn{\overline{w}_j} for this region.
	//' @param x Argument to weight function.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double w_major(const T& x, bool log = true) const {
		static_assert(false, "Specialized implementation for template is needed");
		return 0;
	}

	//' @description
	//' Bifurcate this region into two regions. Use \code{x} as the bifurcation
	//' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
	//' @param x An optional bifurcation point.
	// virtual std::pair<std::unique_ptr<Region<T>>,std::unique_ptr<Region<T>>> bifurcate() const = 0;
	std::pair<Region<T>,Region<T>> bifurcate() const {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}

	// virtual std::pair<std::unique_ptr<Region<T>>,std::unique_ptr<Region<T>>> bifurcate(const T& x) const = 0;
	// virtual std::pair<Region<T>,Region<T>> bifurcate(const T& x) const = 0;
	std::pair<Region<T>,Region<T>> bifurcate(const T& x) const {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}

	Region<T> singleton(const T& x) const {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}

	//' @description
	//' Return a logical value indicating whether this region is bifurcatable.
	bool is_bifurcatable() const {
		static_assert(false, "Specialized implementation for template is needed");
		return false;
	}

	//' @description
	//' The quantity \eqn{\overline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double get_xi_upper(bool log = true) const {
		static_assert(false, "Specialized implementation for template is needed");
		return 0;
	}

	//' @description
	//' The quantity \eqn{\underline{\xi}_j} for this region.
	//' @param log logical; if \code{TRUE}, return result on the log-scale.
	double get_xi_lower(bool log = true) const {
		static_assert(false, "Specialized implementation for template is needed");
		return 0;
	}

	//' @description
	//' A string that describes the region.
	std::string description() const {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}

	//' @description
	//' Print a description of the region.
	void print() const {
		static_assert(false, "Specialized implementation for template is needed");
	};

	/*
	bool operator<(const Region<T>& x) const {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}

	bool operator==(const Region<T>& x) const {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}

	Region<T> operator=(const Region<T>& x) {
		static_assert(false, "Specialized implementation for template is needed");
		return NULL;
	}
	*/
};

}

#endif
