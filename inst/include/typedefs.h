#ifndef VWS_TYPEDEFS_H
#define VWS_TYPEDEFS_H

#include <Rcpp.h>

namespace vws {

static const unsigned int uint_max = std::numeric_limits<unsigned int>::max();

typedef std::function<double(double, bool)> uv_weight_function;
typedef std::function<double(const Rcpp::NumericVector&, bool)> mv_weight_function;

enum class error_action {
	STOP,
	WARNING,
	MESSAGE,
	NONE
};

}

#endif
