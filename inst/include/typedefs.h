#ifndef VWS_TYPEDEFS_H
#define VWS_TYPEDEFS_H

#include <Rcpp.h>
#include "fntl.h"
#include "optimize-hybrid.h"

namespace vws {

/*
* Typedefs for several functions used in VWS programming.
*/
typedef std::function<double(double, bool)> dfdb;
typedef std::function<double(const dfdb& w, double lo, double hi,
	bool log)> optimizer;
typedef std::function<double(double a, double b)> midpoint;

}

#endif
