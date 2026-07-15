// [[Rcpp::depends(RcppArmadillo, vws, fntl)]]
#include <RcppArmadillo.h>
#include <vws.h>
#include <chrono>
#include "ln-norm-proposal.h"
#include "mvnorm.h"
#include "invgamma.h"
#include "logger.h"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs_cpp(const arma::vec& z, const arma::vec& lambda,
	const arma::mat& X, const Rcpp::List& init, const Rcpp::List& prior,
	const Rcpp::List& control, const Rcpp::List& fixed)
{
	unsigned int n = z.n_elem;
	unsigned int d = X.n_cols;

	if (X.n_rows != n) { Rcpp::stop("X.n_rows == n"); }
	if (lambda.n_elem != n) { Rcpp::stop("lambda2.n_elems != n"); }

	const arma::mat& XtX = X.t() * X;

	if (!control.inherits("control_gibbs")) {
		Rcpp::stop("control does not inherit from 'control_gibbs'");
	}
	unsigned int R = control["R"];
	unsigned int burn = control["burn"];
	unsigned int thin = control["thin"];
	unsigned int report = control["report"];
	const arma::uvec& save_latent = control["save_latent"];
	const Rcpp::List& inner_ctrl = control["inner"];

	const Rcpp::String& inner_method = inner_ctrl["method"];
	unsigned int max_rejects = inner_ctrl["max_rejects"];
	double tol_suff = inner_ctrl["tol_suff"];
	double tol_merge = inner_ctrl["tol_merge"];
	unsigned int tune = inner_ctrl["tune"];
	unsigned int N = inner_ctrl["N"];

	unsigned int rep_keep = 0;
	unsigned int R_keep = std::ceil((R - burn) / double(thin));

	// Set up histories
	arma::mat beta_hist(R_keep, d);
	arma::vec sigma2_hist(R_keep);
	arma::mat y_hist(R_keep, save_latent.size());

	// Some statistics on tuning
	arma::uvec rejects_hist(R);
	arma::uvec comps_hist(R);
	arma::uvec tunes_hist(R);
	arma::uvec tuned_hist(R);
	tunes_hist.fill(0);
	tuned_hist.fill(0);
	double avg_comps = 0;

    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	// Set up initial values
	if (!init.inherits("init_gibbs")) {
		Rcpp::stop("init does not inherit from 'init_gibbs'");
	}

	arma::vec beta = init["beta"];
	arma::vec y = init["y"];
	double sigma2 = init["sigma2"];

	arma::vec Xbeta = X * beta;

	// Initialize self-tuned VWS proposals
	std::vector<ln_norm_proposal> proposals;
	for (unsigned int i = 0; i < n; i++) {
	 	ln_norm_proposal x(z(i), lambda(i), Xbeta(i), std::sqrt(sigma2));
	 	proposals.push_back(x);
	}

	// Set up fixed parameters
	if (!fixed.inherits("fixed_gibbs")) {
		Rcpp::stop("fixed does not inherit from 'fixed_gibbs'");
	}

	// Set up prior
	if (!prior.inherits("prior_gibbs")) {
		Rcpp::stop("prior does not inherit from 'prior_gibbs'");
	}

	double sigma_beta = prior["sigma_beta"];
	double a_sigma = prior["a_sigma"];
	double b_sigma = prior["b_sigma"];
	double sigma2_beta = sigma_beta*sigma_beta;

	// Set up timers
	double elapsed_beta = 0;
	double elapsed_sigma2 = 0;
	double elapsed_y = 0;

	for (unsigned int rep = 0; rep < R; rep++)
	{
		// Draw [y | rest]
		if (!fixed["y"]) {
			// Rprintf("Begin Draw [y | rest]\n");
			auto st = std::chrono::system_clock::now();

			if (strcmp(inner_method.get_cstring(), "vws-tune") == 0) {
				// Self-tuned VWS
				for (unsigned int i = 0; i < n; i++) {
					proposals[i].update(Xbeta(i), std::sqrt(sigma2));

					if (rep < tune) {
						const auto& vws_out = vws::rejection_tune(proposals[i], 1, args);
						y(i) = vws_out.draws[0];
						rejects_hist(rep) += vws_out.rejects[0];
						tunes_hist(rep) += vws_out.tunes[0];
						tuned_hist(rep) += (vws_out.tunes[0] > 0);
					} else {
						const auto& vws_out = vws::rejection(proposals[i], 1, args);
						y(i) = vws_out.draws[0];
						rejects_hist(rep) += vws_out.rejects[0];
					}
				}
			} else if (strcmp(inner_method.get_cstring(), "vws-basic") == 0) {
				// VWS without tuning
				for (unsigned int i = 0; i < n; i++) {
					ln_norm_proposal h(z(i), lambda(i), Xbeta(i), std::sqrt(sigma2));
					h.refine(N - 1, tol_suff);
					const auto& vws_out = vws::rejection(h, 1, args);
					y(i) = vws_out.draws[0];
					rejects_hist(rep) += vws_out.rejects[0];
                }
            } else {
                Rcpp::stop("Unrecognized method in inner_ctrl");
            }

			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_y += td.count() * SEC_PER_MICROSEC;
			// Rprintf("End Draw [y | rest]\n");
		}

		// Draw [beta | rest]
		if (!fixed["beta"]) {
			// Rprintf("Begin Draw [beta | rest]\n");
			auto st = std::chrono::system_clock::now();
			const arma::mat& Omega = (1 / sigma2) * XtX + (1 / sigma2_beta) * arma::eye(d,d);
			const arma::vec& mm = arma::solve(Omega, X.t() * arma::log(y) / sigma2);
			beta = r_mvnorm_prec(mm, Omega);
			Xbeta = X * beta;
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_beta += td.count() * SEC_PER_MICROSEC;
			// Rprintf("End Draw [beta | rest]\n");
		}

		// Draw [sigma2 | rest]
		if (!fixed["sigma2"]) {
			// Rprintf("Begin Draw [sigma2 | rest]\n");
			auto st = std::chrono::system_clock::now();
			const arma::vec& e = arma::log(y) - Xbeta;
			double aa = a_sigma + n / 2.0;
			double bb = b_sigma + 1 / 2.0 * dot(e, e);
			sigma2 = r_invgamma(aa, bb);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_sigma2 += td.count() * SEC_PER_MICROSEC;
			// Rprintf("Begin Draw [sigma2 | rest]\n");
		}

		// Save total number of mixture components at this point
		comps_hist(rep) = 0;
		for (unsigned int i = 0; i < n; i++) {
		 	comps_hist(rep) += proposals[i].size();
		}
		avg_comps = comps_hist(rep) / double(n);

		if (rep >= burn && rep % thin == 0) {
			beta_hist.row(rep_keep) = beta.t();
			sigma2_hist[rep_keep] = sigma2;

			for (unsigned int l = 0; l < save_latent.size(); l++) {
				unsigned int i = save_latent(l);
				y_hist(rep_keep, l) = y(i);
			}

			rep_keep++;
		}

		if ((rep + 1) % report == 0) {
			if (strcmp(inner_method.get_cstring(), "vws-tune") == 0) {
				unsigned int s = (rep >= report) ? rep - report : 0;
				unsigned int rejects = arma::sum(rejects_hist(arma::span(s, rep)));
				unsigned int tunes = arma::sum(tunes_hist(arma::span(s, rep)));
				logger("[%d] avg-N: %0.4f  tunes: %d  rejects: %d\n", rep + 1,
					avg_comps, tunes, rejects);
			} else {
				unsigned int s = (rep >= report) ? rep - report : 0;
				unsigned int rejects = arma::sum(rejects_hist(arma::span(s, rep)));
				logger("[%d] rejects: %d\n", rep + 1, rejects);
			}
		}

		Rcpp::checkUserInterrupt();

	}

	Rcpp::List elapsed = Rcpp::List::create(
		Rcpp::Named("beta") = elapsed_beta,
		Rcpp::Named("sigma2") = elapsed_sigma2,
		Rcpp::Named("y") = elapsed_y
	);

	return Rcpp::List::create(
		Rcpp::Named("beta") = beta_hist,
		Rcpp::Named("sigma2") = sigma2_hist,
		Rcpp::Named("y") = y_hist,
		Rcpp::Named("R_keep") = R_keep,
		Rcpp::Named("elapsed") = elapsed,
		Rcpp::Named("R") = R,
		Rcpp::Named("burn") = burn,
		Rcpp::Named("thin") = thin,
		Rcpp::Named("rejects") = rejects_hist,
		Rcpp::Named("comps") = comps_hist,
		Rcpp::Named("tunes") = tunes_hist,
		Rcpp::Named("tuned") = tuned_hist,
		Rcpp::Named("n") = n
	);
}
