library(vws)

source(system.file("examples/lnorm-norm-const-v1/functions.R", package = "vws"))

set.seed(1234)

# ----- Generate data -----
mu = 5
sigma2 = 0.5
lambda2 = 100

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

helper = normal_helper(mean = z, sd = sqrt(lambda2))

support = UnivariateConstRegion$new(a = 0, b = Inf, w = w, g = helper)
regions = list(support)
h_init = FMMProposal$new(regions)

# ----- Refine proposal -----
refine_out = refine(h_init, N = 30, tol = log(1/10), report = 10)
h = refine_out$h

gg = data.frame(bdd = exp(refine_out$log_bdd_hist)) %>%
	mutate(step = row_number() - 1) %>%
	ggplot() +
	geom_line(aes(x = step, y = bdd)) +
	scale_y_continuous(n.breaks = 10) +
	xlab("Step") +
	ylab("Bound") +
	theme_minimal()
print(gg)

bdd = tail(exp(refine_out$log_bdd_hist), 1)
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

# ----- Plot interval based on [y | z] -----
interval_lo = quantile(y, probs = 0.025)
interval_hi = quantile(y, probs = 0.975)
gg + annotate("rect", xmin = interval_lo, xmax = interval_hi, ymin = 0,
		ymax = Inf, alpha = 0.1, fill = "blue") +
	geom_vline(xintercept = z, col = "blue", lwd = 1.05) +
	geom_vline(xintercept = y_true, col = "red", lwd = 1.05)

