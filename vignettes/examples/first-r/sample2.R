sample2 = function(n, kappa, d, N, tol = 0.01, max_rejects = 4*n, report = n / 10)
{
    w = \(x, log = T) {
        out = 0.5 * (d - 3) * log1p(-x^2) + kappa*x
        if (log) { return(out) } else { return(exp(out)) }
    }

    df = \(x, log = F) { dunif(x, min = -1, max = 1, log) }
    pf = \(q, lower.tail = T, log.p = F) {
        punif(q, min = -1, max = 1, lower.tail, log.p)
    }
    qf = \(p, lower.tail = T, log.p = F) {
        qunif(p, min = -1, max = 1, lower.tail, log.p)
    }
    helper = univariate_helper(df, pf, qf)

    maxopt = \(w, lo, hi, log) {
    	A = 1; B = (d-3) / kappa; C = -1
    	x_roots = ( -B + c(-1,1) * sqrt(B^2 - 4*A*C) ) / (2*A)
    	x_root = x_roots[-1 < x_roots & x_roots <= 1]
		if (x_root < lo) {
			out = w(lo, TRUE)
		} else if (x_root > hi) {
			out = w(hi, TRUE)
		} else {
			out = w(x_root, TRUE)
		}

    	return (out)
    }

    minopt = \(w, lo, hi, log) {
		w(min(lo, hi), TRUE)
    }

    ctrl = rejection_control(N = N - 1, tol = tol, max_rejects = max_rejects,
    	report = report, maxopt = maxopt, minopt = minopt)
    rejection_numeric(n = n, lo = -1, hi = 1, w, helper, control = ctrl)
}
