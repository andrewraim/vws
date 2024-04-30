// [[Rcpp::depends(vws)]]
#include "vws.h"

int hello(double x)
{
	printf("Hello world %g\n", x);
	return 0;
}
