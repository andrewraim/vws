library(vws)
library(tidyverse)

source("../common/plots.R")
source("functions.R")
Rcpp::sourceCpp("ln-norm-v1.cpp")
Rcpp::sourceCpp("ln-norm-v2.cpp")
Rcpp::sourceCpp("ln-norm-v3.cpp")

set.seed(1234)

n = 10000

# ----- Generate data -----
mu = 5
sigma = sqrt(0.5)
lambda = 100

y_true = rlnorm(1, mu, sigma)
z = rnorm(1, y_true, lambda)

# ----- Compute the (normalized) target density -----
d_target = get_target_density(mu, sigma, z, lambda)

# List of args for density
args = list(mu = mu, sigma = sigma, z = z, lambda = lambda, log = FALSE)

# ----- Version 1 -----
out = r_ln_norm_v1(n = n, z = z, mu = mu, sigma = sigma,
	lambda = lambda, N = 10, max_rejects = n * 10, report = n / 10)
vws::printf("Empirical rejection rate: %g%%\n",
	100 * sum(out$rejects) / (sum(out$rejects) + n))

plot_density(out$draws) +
	geom_function(fun = d_target, args = args, lty = 2)
plot_bounds(out$lbdd)

# ----- Version 2 -----
out = r_ln_norm_v2(n = n, z = z, mu = mu, sigma = sigma,
	lambda = lambda, N = 10, max_rejects = n * 10, report = n / 10)
vws::printf("Empirical rejection rate: %g%%\n",
	100 * sum(out$rejects) / (sum(out$rejects) + n))

plot_density(out$draws) +
	geom_function(fun = d_target, args = args, lty = 2)
plot_bounds(out$lbdd)

# ----- Version 3 -----
out = r_ln_norm_v3(n = n, z = z, mu = mu, sigma = sigma,
	lambda = lambda, N = 10, max_rejects = n * 10, report = n / 10)
vws::printf("Empirical rejection rate: %g%%\n",
	100 * sum(out$rejects) / (sum(out$rejects) + n))

plot_density(out$draws) +
	geom_function(fun = d_target, args = args, lty = 2)
plot_bounds(out$lbdd)
