library(vws)
library(tidyverse)

source("functions.R")
Rcpp::sourceCpp("lnorm-norm-v1.cpp")

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

# ----- Version 1 -----
out = r_ln_norm_v1(n = n, z = z, mu = mu, sigma = sigma,
	lambda = lambda, N = 10, max_rejects = n * 10, report = n / 10)
vws::printf("Empirical rejection rate: %g%%\n",
	100 * sum(out$rejects) / (sum(out$rejects) + n))

gg = ggplot() +
	geom_function(fun = d_target, args = list(log = FALSE)) +
	geom_density(data = data.frame(x = out$draws), aes(x = x), col = "blue", lty = 2) +
	xlab("x") +
	ylab("Density") +
	theme_linedraw()
print(gg)

# ----- Version 2 -----
out = r_ln_norm_v2(n = n, z = z, mu = mu, sigma = sigma,
	lambda = lambda, N = 10, max_rejects = n * 10, report = n / 10)
vws::printf("Empirical rejection rate: %g%%\n",
	100 * sum(out$rejects) / (sum(out$rejects) + n))

gg = ggplot() +
	geom_function(fun = d_target, args = list(log = FALSE)) +
	geom_density(data = data.frame(x = out$draws), aes(x = x), col = "blue", lty = 2) +
	xlab("x") +
	ylab("Density") +
	theme_linedraw()
print(gg)

