library(vws)

source(system.file("examples/lnorm-norm-const-v1/functions.R", package = "vws"))
source(system.file("examples/lnorm-norm-linear/CustomLinearRegion.R", package = "vws"))

# ----- Generate data -----
mu = 5
sigma2 = 0.5
lambda2 = 100

set.seed(1234)
y_true = rlnorm(1, mu, sqrt(sigma2))
z = rnorm(1, y_true, sqrt(lambda2))
print(y_true)
print(z)

# ----- Construct proposal -----
w = function(y, log = TRUE) {
	out = -log(y) - (log(y) - mu)^2 / (2*sigma2) + log(y > 0)
	out[y == 0] = -Inf
	if (log) { return(out) } else { return(exp(out)) }
}

y_star = exp(mu - sigma2 + 1)
region1 = CustomLinearRegion$new(a = 1e-6, b = y_star, mu, sigma2, z, lambda2)
region2 = CustomLinearRegion$new(a = y_star, b = 1e6, mu, sigma2, z, lambda2)
regions = list(region1, region2)

h_init = FMMProposal$new(regions)

# ----- Adapt proposal -----
adapt_out = adapt(h_init, N = 30, report = 10)
h = adapt_out$h

gg = data.frame(bdd = exp(adapt_out$log_bdd_hist)) %>%
	mutate(step = row_number() - 1) %>%
	ggplot() +
	geom_line(aes(x = step, y = bdd)) +
	scale_y_continuous(n.breaks = 10) +
	xlab("Step") +
	ylab("Bound") +
	theme_minimal()
print(gg)

bdd = tail(exp(adapt_out$log_bdd_hist), 1)
cat("Upper bound for percent of rejections:", 100 * bdd)

# ----- Rejection sampling -----
ctrl = rejection_control(report = 5000, extra_outputs = TRUE)
out = rejection(h, n = 10000, control = ctrl)
y = unlist(out$draws)

cat("Percent of proposed draws which were rejected:",
	sum(out$rejects) / (length(y) + sum(out$rejects)) * 100)

# ----- Plot draws vs. target density -----
gg = data.frame(y = y) %>%
	ggplot() +
	geom_histogram(aes(x = y, y = after_stat(density)), col = "black",
		fill = NA, bins = 25) +
	geom_function(fun = d_target, args = list(mu = mu, sigma2 = sigma2, z = z,
		lambda2 = lambda2, log = FALSE), lty = 2) +
	xlab("y") +
	ylab("Density") +
	theme_minimal()
print(gg)

# ----- Plot log f(x) vs. log h(x) -----
log_h = function(y) { sapply(y, h$d, log = TRUE) }
log_f = function(y) { d_target(y, mu, sigma2, z, lambda2, log = TRUE) }

xlim = range(y)
ggplot() +
	geom_function(fun = log_f, col = "orange", lty = 1) +
	geom_function(fun = log_h, col = "black", lty = 2) +
	scale_x_continuous(limits = xlim) +
	xlab("y") +
	ylab("Density") +
	theme_minimal()
