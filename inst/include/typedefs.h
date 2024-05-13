#ifndef VWS_TYPEDEFS_H
#define VWS_TYPEDEFS_H

#include <Rcpp.h>

namespace vws {

typedef std::function<double(double, bool)> weight_function;
enum class ErrorAction { STOP, WARNING, MESSAGE, NONE };

}

#endif
