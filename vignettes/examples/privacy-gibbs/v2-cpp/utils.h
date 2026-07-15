#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

inline void logger(const char* fmt, ...)
{
	// Get the current time in the local time zone
	std::time_t raw = std::time(nullptr);
	std::tm* local = std::localtime(&raw);

	// Write formatted time to a string
	char buffer[64];
	strftime(buffer, 64, "%Y-%m-%d %H:%M:%S", local);

	// Insert "..." placeholders into format string for message; see
	// <https://stackoverflow.com/q/1056411>
	char msg[256];
	va_list args;
	va_start(args, fmt);
	vsnprintf(msg, 255, fmt, args);
	va_end(args);

	// Print the formatted message with timestamp
	Rprintf("%s - %s", buffer, msg);
}

inline double r_invgamma(double a, double b)
{
    return 1 / R::rgamma(a, 1 / b);
}

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
