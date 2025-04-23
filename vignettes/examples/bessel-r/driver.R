library(tidyverse)
library(vws)

source("bessel.R")
source("sample.R")
source("sample2.R")

n = 20000
a = 5
nu = 2
N = 12
tol = 0
max_rejects = 50000
report = 10000

# ----- Version 1 -----
# Use numerical optimization to compute constants in majorizer
out = sample(n, a, nu, N, tol, max_rejects, report)

xseq = seq(0, max(out$draws))
fseq = d_bessel(xseq, a = 5, nu = 2)

data.frame(x = out$draws) %>%
	ggplot() +
	geom_bar(aes(x, y = after_stat(count / sum(count))), fill = "white", col = "black") +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y)) +
	xlab("x") +
	ylab("Empirical Density") +
	theme_light()

data.frame(lbdd = out$lbdd) %>%
	mutate(step = row_number()) %>%
	ggplot() +
	geom_line(aes(step, lbdd)) +
	theme_light()

# ----- Version 2 -----
# Use custom optimization routine to compute constants in majorizer
out2 = sample2(n, a, nu, N, tol, max_rejects, report)

xseq = seq(0, max(out2$draws))
fseq = d_bessel(xseq, a, nu)

data.frame(x = out2$draws) %>%
	ggplot() +
	geom_bar(aes(x, y = after_stat(count / sum(count))), fill = "white", col = "black") +
	geom_point(data = data.frame(x = xseq, y = fseq), aes(x,y)) +
	xlab("x") +
	ylab("Empirical Density") +
	theme_light()

data.frame(lbdd = out2$lbdd) %>%
	mutate(step = row_number()) %>%
	ggplot() +
	geom_line(aes(step, lbdd)) +
	theme_light()

