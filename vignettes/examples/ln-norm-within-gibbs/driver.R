library(tidyverse)
source("gibbs.R")

# ----- Generate data from the model -----
set.seed(1234)

n = 400
d = 2
x = rnorm(n)
X = cbind(1, x)
sigma_true = 0.5
beta_true = c(1, 2)
Xbeta_true = X %*% beta_true
y_true = rlnorm(n, Xbeta_true, sigma_true)
lambda = rep(0.5, n)
z = rnorm(n, y_true, lambda)

# ----- Fit the model using self-tuned VWS -----
fixed = fixed_gibbs()
init = init_gibbs(n, d, y = y_true)
inner = control_inner(tol_suff = 0.85, tol_merge = 0.01, tune = 100, , method = "vws-tune")
prior = prior_gibbs()
control = control_gibbs(R = 3000, burn = 1000, report = 100, inner = inner,
	save_latent = 1:10)

gibbs_out = gibbs(z, lambda, X, init, prior, control, fixed)
print(gibbs_out)

# ----- Fit the model using basic VWS without tuning -----
# Note this takes longer

fixed = fixed_gibbs(y = FALSE)
init = init_gibbs(n, d, y = y_true)
inner = control_inner(tol_suff = 0.85, method = "vws-basic")
prior = prior_gibbs()
control = control_gibbs(R = 3000, burn = 1000, report = 100, inner = inner,
	save_latent = 1:10)

gibbs0_out = gibbs(z, lambda, X, init, prior, control, fixed)
print(gibbs0_out)

# ----- Plot draws of parameters -----
df_draws = with(gibbs_out,
		data.frame(beta1 = beta[,1], beta2 = beta[,2], sigma2 = sigma2)) %>%
	mutate(iter = row_number())

df_draws %>%
	ggplot() +
	geom_line(aes(x = iter, y = beta1)) +
	xlab("Iteration") +
	ylab(expression(beta[1])) +
	theme_minimal()

df_draws %>%
	ggplot() +
	geom_line(aes(x = iter, y = beta2)) +
	xlab("Iteration") +
	ylab(expression(beta[2])) +
	theme_minimal()

df_draws %>%
	ggplot() +
	geom_line(aes(x = iter, y = sigma2)) +
	xlab("Iteration") +
	ylab(expression(sigma^2)) +
	theme_minimal()

# ----- Plot draws of a selected latent value -----
df_latent = as.data.frame(gibbs_out$y) %>%
	rename_with(~ gsub("V", "y", .x)) %>%
	mutate(iter = row_number())

df_latent %>%
	ggplot() +
	geom_line(aes(x = iter, y = y3)) +
	geom_hline(yintercept = y_true[3], col = "red", lty = 2) +
	xlab("Iteration") +
	ylab(expression(y[3])) +
	theme_minimal()

# ----- Plot tuning metrics per iteration -----
df_metrics = with(gibbs_out, data.frame(rejects, comps, tunes, tuned)) %>%
	mutate(iter = row_number())

df_metrics %>%
	filter(iter > 1) %>%
	ggplot() +
	geom_line(aes(x = iter, y = rejects)) +
	xlab("Iteration") +
	ylab("Rejects") +
	geom_vline(xintercept = inner$tune, col = "red", lty = 2) +
	theme_minimal()

df_metrics %>%
	filter(iter > 1) %>%
	ggplot() +
	geom_line(aes(x = iter, y = tunes)) +
	xlab("Iteration") +
	ylab("Tuning Adjustments") +
	geom_vline(xintercept = inner$tune, col = "red", lty = 2) +
	theme_minimal()
