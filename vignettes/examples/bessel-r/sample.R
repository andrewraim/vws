sample = function(n, a, nu, N, tol = 0.01, max_rejects = 4*n, report = n / 10)
{
   lambda = a^2 / 4

    w = \(x, log = T) {
        out = -lgamma(x + nu + 1)
        ifelse(log, out, exp(out))
    }

    df = \(x, log = F) { dpois(x, lambda, log) }
    pf = \(q, lower.tail = T, log.p = F) { ppois(q, lambda, lower.tail, log.p) }
    qf = \(p, lower.tail = T, log.p = F) { qpois(p, lambda, lower.tail, log.p) }
    sf = \(x) { x >= 0 }
    helper = univariate_helper(df, pf, qf, sf)

    ctrl = rejection_control(N = N - 1, tol = tol, max_rejects = max_rejects,
    	report = report)
    rejection_int(n = n, lo = -Inf, hi = Inf, w, helper, control = ctrl)
}
