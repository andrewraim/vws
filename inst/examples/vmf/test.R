# ----- Test functions for base distribution -----
# Compare empirical draws to functions.

kappa = 0.2
x = r_base(100000, kappa)
df = data.frame(x = x)

x_seq = seq(-1, 1, length.out = 100)
q_seq = seq(0, 1, length.out = 100)

d_out = d_base(x_seq, kappa)
p_out = p_base(x_seq, kappa)
q_out = q_base(q_seq, kappa)
q_emp = quantile(x, probs = q_seq)

# Compare histogram to density function
ggplot(df) +
	geom_histogram(aes(x, after_stat(density)), col = "black",
		fill = "yellow", bins = 25) +
	geom_line(data = data.frame(x = x_seq, y = d_out), aes(x, y)) +
	scale_x_continuous(limits = c(1e-4 - 1, 1 - 1e-4)) +
	theme_minimal()

# Compare ECDF to CDF function
ggplot(df) +
	stat_ecdf(aes(x), geom = "step") +
	geom_function(fun = p_base, args = list(kappa = kappa), col = "blue", lty = 2) +
	theme_minimal()

# Compare sample quantiles to quantile function
df = data.frame(q = q_seq, q_emp = q_emp, q_out = q_out)
ggplot(df) +
	geom_line(aes(q, q_emp)) +
	geom_line(aes(q, q_out), lty = 2, col = "blue") +
	xlab("Probability") +
	xlab("Quantile") +
	theme_minimal()

# ----- Test functions for target distribution -----
kappa = 0.2
d = 2

x_seq = seq(-1, 1, length.out = 100)

nc_target(kappa, d, log = FALSE)
d_seq = d_target(x_seq, kappa, d)
p_seq = p_target(x_seq, kappa, d)

df = data.frame(x = x_seq, d = d_seq, p = p_seq)

ggplot(df) +
	geom_line(aes(x, d)) +
	scale_x_continuous(limits = c(1e-4 - 1, 1 - 1e-4)) +
	xlab("x") +
	ylab("Density") +
	theme_minimal()

ggplot(df) +
	geom_line(aes(x, p)) +
	scale_x_continuous(limits = c(1e-4 - 1, 1 - 1e-4)) +
	xlab("x") +
	ylab("CDF") +
	theme_minimal()
