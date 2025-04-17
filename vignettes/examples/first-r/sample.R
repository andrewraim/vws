library(vws)

sample = function(n, kappa, d, N)
{
	w = \(x, log = T) {
		out = -Inf;
		if (-1 < x && x <= 1) {
			out = 0.5 * (d - 3) * log1p(-x^2) + kappa*x
		}
		if (log) { return(out) } else { return(exp(out)) }
	}

	df = \(x, log = F) {
		dunif(x, min = -1, max = 1, log = log)
	}
	pf = \(q, lower.tail = T, log.p = F) {
		punif(q, min = -1, max = 1, lower.tail = lower.tail, log.p = log)
	}
	qf = \(p, lower.tail = T, log.p = F) {
		qunif(p, min = -1, max = 1, lower.tail = lower.tail, log.p = log)
	}
	sf = \(x) {
		-1 <= x && x <= 1
	}
	helper = univariate_helper(df, pf, qf, sf)

	ctrl = rejection_control(max_rejects = 1000, report = 100)
	rejection(n = n, lo = -1, hi = 1, w, helper, N = N - 1,
		control = ctrl)
}
