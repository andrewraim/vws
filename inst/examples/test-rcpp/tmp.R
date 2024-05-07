library(vws)
Rcpp::sourceCpp("sampler.cpp", rebuild = TRUE)

# ----- Generate data -----
z = 1
mu = 5
sigma2 = 0.5
lambda2 = 100

set.seed(1234)
y_true = rlnorm(1, mu, sqrt(sigma2))
z = rnorm(1, y_true, sqrt(lambda2))

# ----- Sampler in R -----

w = function(y, log = TRUE) {
	out = -log(y) - (log(y) - mu)^2 / (2*sigma2) + log(y > 0)
	out[y == 0] = -Inf
	if (log) { return(out) } else { return(exp(out)) }
}

helper = normal_univariate_helper(mean = z, sd = sqrt(lambda2))

support = UnivariateConstRegion$new(a = 0, b = Inf, w = w, g = helper)
regions = list(support)
h_init = FMMProposal$new(regions)

# adapt_out = adapt(h_init, N = 9, report = 10)
# h = adapt_out$h
# print(h, n = 10)
h = h_init

# bdd = tail(exp(adapt_out$log_bdd_hist), 1)

ctrl = rejection_control(report = 5000, extra_outputs = TRUE)
out = rejection(h, n = 10000, control = ctrl)
y = unlist(out$draws)

# ----- Sampler in C++ -----
r_lognormal_normal(n = 10, z = z, mu = mu, sigma2 = sigma2, lambda2 = lambda2)
