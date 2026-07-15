#ifndef INVGAMMA_H
#define INVGAMMA_H

#include <Rcpp.h>

inline double r_invgamma(double a, double b)
{
    return 1 / R::rgamma(a, 1 / b);
}

#endif
