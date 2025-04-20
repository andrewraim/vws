sample = function(n, kappa, d, N, tol = 0.01, max_rejects = 4*n, report = n / 10)
{
    w = \(x, log = T) {
        out = -Inf
        if (-1 < x && x <= 1) { out = 0.5 * (d - 3) * log1p(-x^2) + kappa*x }
        if (log) { return(out) } else { return(exp(out)) }
    }

    df = \(x, log = F) { dunif(x, min = -1, max = 1, log) }
    pf = \(q, lower.tail = T, log.p = F) {
        punif(q, min = -1, max = 1, lower.tail, log.p)
    }
    qf = \(p, lower.tail = T, log.p = F) {
        qunif(p, min = -1, max = 1, lower.tail, log.p)
    }
    sf = \(x) { -1 <= x & x <= 1 }
    helper = univariate_helper(df, pf, qf, sf)

    ctrl = rejection_control(N = N - 1, tol = tol, max_rejects = max_rejects,
    	report = report)
    rejection_numeric(n = n, lo = -1, hi = 1, w, helper, control = ctrl)
}
