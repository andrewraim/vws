#ifndef MGF_TRUNCNORM_H
#define MGF_TRUNCNORM_H

#include <Rcpp.h>

// MGF needed for choosing expansion point in LinearVWSRegion
// [[Rcpp::export]]
double mgf_truncnorm(double s, double a, double b, double z, double lambda,
	bool log = false);

#endif
