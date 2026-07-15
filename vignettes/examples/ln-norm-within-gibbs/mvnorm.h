#ifndef MVNORM_H
#define MVNORM_H

#include <RcppArmadillo.h>

inline arma::mat r_mvnorm_prec(const arma::vec& mu, const arma::mat& Omega)
{
	if (mu.n_elem != Omega.n_rows) {
		Rcpp::stop("mu.n_elem != Omega.n_rows");
	}

	if (Omega.n_rows != Omega.n_cols) {
		Rcpp::stop("Omega.n_rows != Omega.n_cols");
	}

    unsigned int k = mu.n_elem;
    const arma::vec& z = arma::randn(k);
    const arma::mat& A = arma::chol(Omega);
    return arma::solve(A, z) + mu;
}

#endif
