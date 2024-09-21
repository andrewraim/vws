# Operations for Base Distribution

d_base = function(x, kappa, log = FALSE)
{
	if (kappa < 0) {
		out = -kappa*x + log(kappa) - vws::log_sub2_exp(-kappa,kappa) + log(-1 < x & x < 1)
	} else {
		out = kappa*x + log(kappa) - vws::log_sub2_exp(kappa,-kappa) + log(-1 < x & x < 1)
	}
	if (log) { return(out) } else { return(exp(out)) }
}

p_base = function(x, kappa, lower.tail = TRUE, log.p = FALSE)
{
	n = length(x)
	islower = rep(lower.tail, n)
	lo = ifelse(islower, rep(-1, n), x)
	hi = ifelse(!islower, rep(1, n), x)
	ones = rep(1, n)
	if (kappa < 0) {
		log_p_lo = log_sub2_exp(-kappa*ones, kappa*lo) - log_sub2_exp(-kappa*ones, kappa*ones)
		log_p_hi = log_sub2_exp(-kappa*ones, kappa*hi) - log_sub2_exp(-kappa*ones, kappa*ones)
	} else {
		log_p_lo = log_sub2_exp(kappa*lo, -kappa*ones) - log_sub2_exp(kappa*ones, -kappa*ones)
		log_p_hi = log_sub2_exp(kappa*hi, -kappa*ones) - log_sub2_exp(kappa*ones, -kappa*ones)
	}
	out = log_sub2_exp(log_p_hi, log_p_lo)
	if (log.p) { return(out) } else { return(exp(out)) }
}

q_base = function(p, kappa, lower.tail = TRUE, log.p = FALSE)
{
	n = length(p)
	if (log.p) { lp = p } else { lp = log(p) }
	if (!lower.tail) { lp = log_sub2_exp(rep(0, n), lp) }
	ones = rep(1, n)

	if (kappa < 0) {
		out = 1/kappa * vws::log_sub2_exp(-kappa*ones, lp +
			vws::log_sub2_exp(-kappa*ones, kappa*ones))
	} else {
		out = 1/kappa * vws::log_add2_exp(-kappa*ones, lp +
			vws::log_sub2_exp(kappa*ones, -kappa*ones))
	}
	return(out)
}

r_base = function(n, kappa)
{
	u = runif(n)
	q_base(u, kappa)
}

vmf_base_helper = function(kappa)
{
	univariate_helper(
		d = function(x, log = FALSE) {
			d_base(x, kappa = kappa, log = log)
		},
		p = function(q, lower.tail = TRUE, log.p = FALSE) {
			p_base(q, kappa = kappa, lower.tail = lower.tail, log.p = log.p)
		},
		q = function(p, lower.tail = TRUE, log.p = FALSE) {
			q_base(p, kappa = kappa, lower.tail = lower.tail, log.p = log.p)
		},
		s = function(x) { -1 < x & x < 1 }
	)
}
