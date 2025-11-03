#ifndef VWS_H
#define VWS_H

/*
* Export classes, functions, etc to package users.
*
* Note that `result.h` must be listed first because it uses a technique
* described in [1] to serialize objects between R and C++.
*
* [1] Dirk Eddelbuettel and Romain Francois. Rcpp Extending. Vignette from
* Rcpp version 1.0.13.
*/

#include "result.h"

#include "categ.h"
#include "FMMProposal.h"
#include "gumbel.h"
#include "IntConstRegion.h"
#include "log-sum-exp.h"
#include "logit.h"
#include "optimize-hybrid.h"
#include "Region.h"
#include "RealConstRegion.h"
#include "RealConstRegion-defaults.h"
#include "rect.h"
#include "rejection.h"
#include "rejection-args.h"
#include "seq.h"
#include "timestamp.h"
#include "typedefs.h"
#include "UnivariateHelper.h"

#endif

