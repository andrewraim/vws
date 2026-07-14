#include <RcppArmadillo.h>
#include <chrono>
#include "saevws.h"
#include "armspp"

const double SEC_PER_MICROSEC = 1e-6;

// [[Rcpp::export]]
Rcpp::List gibbs(const arma::vec& y, const arma::vec& lambda,
	const arma::mat& X, const Rcpp::List& init, const Rcpp::List& control,
	const Rcpp::List& fixed)
{
	unsigned int m = y.n_elem;
	unsigned int d = X.n_cols;

	stopifnot(X.n_rows == m, "X.n_rows == m");
	stopifnot(lambda.n_elem == m, "lambda2.n_elems == m");

	const arma::mat& XtX = crossprod(X);

	stopifnot(control.inherits("control"), "control inherits from control");
	unsigned int R = control["R"];
	unsigned int burn = control["burn"];
	unsigned int thin = control["thin"];
	unsigned int report = control["report"];
	const arma::uvec& save_latent = control["save_latent"];
	const Rcpp::List& inner_ctrl = control["inner"];

	unsigned int max_rejects = inner_ctrl["max_rejects"];
	double tol_suff = inner_ctrl["tol_suff"];
	double tol_merge = inner_ctrl["tol_merge"];
	unsigned int N = inner_ctrl["N"];
	unsigned int tune = inner_ctrl["tune"];

	unsigned int rep_keep = 0;
	unsigned int R_keep = std::ceil((R - burn) / double(thin));

	// Set up histories
	arma::mat beta_hist(R_keep, d);
	arma::vec sigma2_hist(R_keep);
	arma::mat mu_hist(R_keep, save_latent.size());

	arma::uvec rejects_hist(R);
	arma::uvec comps_hist(R);
	arma::uvec tunes_hist(R);
	arma::uvec tuned_hist(R);
	arma::uvec rejects_areas(m);
	tunes_hist.fill(0);
	tuned_hist.fill(0);
	rejects_areas.fill(0);
	double avg_comps = 0;

    vws::rejection_args args;
    args.max_rejects = max_rejects;
    args.report = 1e6;
	args.tol_suff = tol_suff;
	args.tol_merge = tol_merge;

	// Set up initial values
	stopifnot(init.inherits("init"), "init inherits from init");
	arma::vec beta = init["beta"];
	arma::vec mu = init["mu"];
	double sigma2 = init["sigma2"];

	arma::vec Xbeta = X * beta;

	// Initialize self-tuned VWS proposals
	std::vector<ln_norm_proposal> proposals;
	for (unsigned int i = 0; i < m; i++) {
	 	ln_norm_proposal x(y(i), lambda(i), Xbeta(i), std::sqrt(sigma2));
	 	proposals.push_back(x);
	}

	// Set up fixed parameters
	stopifnot(fixed.inherits("fixed"), "fixed inherits from fixed");

	// Note: flat prior is assumed for each parameter in this model

	// Set up timers
	double elapsed_beta = 0;
	double elapsed_sigma2 = 0;
	double elapsed_mu = 0;

	for (unsigned int rep = 0; rep < R; rep++)
	{
		// Draw [mu | rest]
		if (!fixed["mu"]) {
			auto st = std::chrono::system_clock::now();

			/*
			* Self-tuned VWS using vws package
			*/
			for (unsigned int i = 0; i < m; i++) {
				proposals[i].update(Xbeta(i), std::sqrt(sigma2));

				if (rep < tune) {
					const auto& vws_out = vws::rejection_tune(proposals[i], 1, args);
					mu(i) = vws_out.draws[0];
					rejects_areas(i) += vws_out.rejects[0];
					rejects_hist(rep) += vws_out.rejects[0];
					tunes_hist(rep) += vws_out.tunes[0];
					tuned_hist(rep) += (vws_out.tunes[0] > 0);
				} else {
					const auto& vws_out = vws::rejection(proposals[i], 1, args);
					mu(i) = vws_out.draws[0];
					rejects_areas(i) += vws_out.rejects[0];
					rejects_hist(rep) += vws_out.rejects[0];
				}
			}

			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_mu += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [beta | rest]
		if (!fixed["beta"]) {
			auto st = std::chrono::system_clock::now();
			const arma::vec& mm = arma::solve(XtX, crossprod(X, arma::log(mu)));
			const arma::mat& Omega = (1 / sigma2) * XtX;
			beta = r_mvnorm_prec(mm, Omega);
			Xbeta = X * beta;
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_beta += td.count() * SEC_PER_MICROSEC;
		}

		// Draw [sigma2 | rest]
		if (!fixed["sigma2"]) {
			auto st = std::chrono::system_clock::now();
			double aa = m / 2.0;
			double bb = 1 / 2.0 * dot(log(mu) - Xbeta);
			sigma2 = r_invgamma(aa, bb);
			auto et = std::chrono::system_clock::now();
			auto td = std::chrono::duration_cast<std::chrono::microseconds>(et - st);
			elapsed_sigma2 += td.count() * SEC_PER_MICROSEC;
		}

		// Save total number of mixture components at this point
		comps_hist(rep) = 0;
		for (unsigned int i = 0; i < m; i++) {
		 	comps_hist(rep) += proposals[i].size();
		}
		avg_comps = comps_hist(rep) / double(m);

		if (rep >= burn && rep % thin == 0) {
			beta_hist.row(rep_keep) = beta.t();
			sigma2_hist[rep_keep] = sigma2;

			for (unsigned int l = 0; l < save_latent.size(); l++) {
				unsigned int i = save_latent(l);
				mu_hist(rep_keep, l) = mu(i);
			}

			rep_keep++;
		}

		if ((rep + 1) % report == 0) {
			unsigned int s = (rep >= report) ? rep - report : 0;
			unsigned int rejects = arma::sum(rejects_hist(arma::span(s, rep)));
			unsigned int tunes = arma::sum(tunes_hist(arma::span(s, rep)));
			logger("[%d] avg-N: %0.4f  tunes: %d  rejects: %d\n", rep + 1,
				avg_comps, tunes, rejects);
		}

		Rcpp::checkUserInterrupt();

	}

	Rcpp::List elapsed = Rcpp::List::create(
		Rcpp::Named("beta") = elapsed_beta,
		Rcpp::Named("sigma2") = elapsed_sigma2,
		Rcpp::Named("mu") = elapsed_mu
	);

	return Rcpp::List::create(
		Rcpp::Named("beta") = beta_hist,
		Rcpp::Named("sigma2") = sigma2_hist,
		Rcpp::Named("mu") = mu_hist,
		Rcpp::Named("R_keep") = R_keep,
		Rcpp::Named("elapsed") = elapsed,
		Rcpp::Named("R") = R,
		Rcpp::Named("burn") = burn,
		Rcpp::Named("thin") = thin,
		Rcpp::Named("inner_method") = inner_method,
		Rcpp::Named("rejects") = rejects_hist,
		Rcpp::Named("rejects_areas") = rejects_areas,
		Rcpp::Named("comps") = comps_hist,
		Rcpp::Named("tunes") = tunes_hist,
		Rcpp::Named("tuned") = tuned_hist,
		Rcpp::Named("m") = m,
		Rcpp::Named("inner_method") = inner_method
	);
}
