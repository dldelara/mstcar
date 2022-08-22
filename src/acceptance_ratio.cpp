#include <RcppArmadillo.h>
#include "save_load.h"
#include "update_helper.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::cube acpt_cube(arma::field<arma::cube> f, int Ng, int Nt, int Ns) {
	int l = f.n_elem;
	cube diff(Ng, Nt, Ns, fill::zeros);
	cube c;
	for (int i = 0; i < l - 2; i++) {
		c = f[i] - f[i + 1];
		for (unsigned int j = 0; j < c.n_elem; j++) if (c[j] != 0) diff[j]++;
	}
	diff /= (l - 1);
	return diff;
}
//[[Rcpp::export]]
arma::rowvec acpt_vec(arma::field<arma::mat> f, int Ng) {
	int l = f.n_elem;
	rowvec diff(Ng, fill::zeros);
	rowvec m;
	for (int i = 0; i < l - 2; i++) {
		m = f[i] - f[i + 1];
		for (unsigned int j = 0; j < m.n_elem; j++) if (m[j] != 0) diff[j]++;
	}
	diff /= (l - 1);
	return diff;
}
//[[Rcpp::export]]
arma::cube acceptance_ratio_cube(Rcpp::List mod, arma::vec file_suff, int burn) {
	List params = mod["params"];
	String name = params["name"];
	String dir  = params["dir"];
	vec dNd     = params["dNd"];
	int Ng      = dNd[0];
	int Nt      = dNd[1];
	int Ns      = dNd[2];
	List sample = get_output(mod, burn, 1, "theta", file_suff);
	field<cube> theta = sample["theta"];
	return acpt_cube(theta, Ng, Nt, Ns);
}
//[[Rcpp::export]]
arma::rowvec acceptance_ratio_mat(Rcpp::List mod, arma::vec file_suff, int burn) {
	List params = mod["params"];
	String name = params["name"];
	String dir  = params["dir"];
	vec dNd     = params["dNd"];
	int Ng      = dNd[0];
	List sample = get_output(mod, burn, 1, "rho", file_suff);
	field<mat> rho = sample["rho"];
	return acpt_vec(rho, Ng);
}
//[[Rcpp::export]]
void adaptive_variance_theta(arma::cube& gamma, arma::field<arma::cube> new_theta, int s, int Ng, int Nt, int Ns) {
	field<cube> theta_a = new_theta.rows(s - 99, s);
	cube gam_acpt = acpt_cube(theta_a, Ng, Nt, Ns);
	for (int i = 0; i < Ng * Nt * Ns; i++) {
		if (gam_acpt[i] > 0.75) gam_acpt[i] = 0.75;
		if (gam_acpt[i] < 0.2) gam_acpt[i] = 0.2;
		gamma[i] *= gam_acpt[i] / 0.43;
	}
}
//[[Rcpp::export]]
void adaptive_variance_rho(arma::rowvec& delta, arma::field<arma::mat> new_rho, int s, int Ng) {
	field<mat> rho_a = new_rho.rows(s - 99, s);
	rowvec rho_acpt = acpt_vec(rho_a, Ng);
	for (int k = 0; k < Ng; k++) {
		if (rho_acpt[k] > 0.75) rho_acpt[k] = 0.75;
		if (rho_acpt[k] < 0.2) rho_acpt[k] = 0.2;
		delta[k] *= rho_acpt[k] / 0.43;
	}
}