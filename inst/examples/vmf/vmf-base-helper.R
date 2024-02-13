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
