library(vws)
library(statmod)

d_target = function(x, mu, sigma, z, lambda, log = FALSE, normalize = TRUE, Q = 10)
{
	psi = 1
	if (normalize) {
		## Compute normalizing constant using Hermite quadrature with Q points
		w = function(x, log = TRUE) {
			out = -log(x) - (log(x) - mu)^2 / (2*sigma^2) + log(x > 0)
			out[x == 0] = -Inf
			if (log) { return(out) } else { return(exp(out)) }
		}

		q = function(x, log = FALSE) {
			tx = z - sqrt(2) * lambda * x
			ii = which(tx > 0)
			out = rep(-Inf, length(x))
			out[ii] = -log(tx[ii]) - (log(tx[ii]) - mu)^2 / (2*sigma^2) - 1/2 * log(pi)
			if (log) { return(out) } else { return(exp(out)) }
		}

		quad_out = gauss.quad(n = Q, kind = "hermite")
		ww = quad_out$weights
		xx = quad_out$nodes
		psi = sum(ww * q(xx))
	}

	out = w(x, log = TRUE) + dnorm(x, z, lambda, log = TRUE) - log(psi)
	if (log) { return(out) } else { return(exp(out)) }
}

plot_weight_function = function(mu, sigma)
{
	w = function(y, log = TRUE) {
		out = -log(y) - (log(y) - mu)^2 / (2*sigma^2) + log(y > 0)
		out[y == 0] = -Inf
		if (log) { return(out) } else { return(exp(out)) }
	}

	xlim = c(30, 280)
	yseq = seq(xlim[1], xlim[2], length.out = 200)
	wy = w(yseq, log = FALSE)
	wmin = min(wy)

	y1 = exp(mu - sigma^2)
	y2 = exp(mu - sigma^2 + 1)

	w1 = w(y1, log = FALSE)
	w2 = w(y2, log = FALSE)

	ggplot(data.frame(y = yseq, wy = wy)) +
		geom_line(aes(y, wy)) +
		geom_vline(xintercept = y1, lty = 2, col = "blue", lwd = 1.0) +
		geom_vline(xintercept = y2, lty = 2, col = "brown", lwd = 1.0) +
		scale_x_continuous(expand = c(0,0), n.breaks = 10) +
		xlab("y") +
		ylab("w(y)") +
		theme_light() +
		theme(
			axis.line.x = element_line(arrow = grid::arrow(
				length = unit(0.2, "cm"), ends = "last", type = "closed")),
			axis.line.y = element_line(arrow = grid::arrow(
				length = unit(0.2, "cm"), ends = "last", type = "closed"))
	)
}
