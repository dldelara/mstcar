#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "update_helper.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_theta(arma::cube& theta,
	arma::cube Y, arma::cube n,                      // data
	arma::rowvec tau2, arma::mat beta, arma::cube z, // parameters
	arma::cube gamma,                                // hyperparameters
	int Ng, int Nt, int Ns, Rcpp::String method      // metaparameters
) {
	for (int i = 0; i < Ng * Nt * Ns; i++) {
		double theta_star = R::rnorm(theta[i], gamma[i]);
		double r_t = Y[i] * (theta_star - theta[i]);
		if (method == "binom") r_t += - n[i] * (log(1 + exp(theta_star)) - log(1 + exp(theta[i])));
		if (method == "pois" ) r_t += - n[i] * (exp(theta_star) - exp(theta[i]));
		r_t += - 1 / (2 * tau2[i % Ng]) * (pow(theta_star - beta[i % (Ng * Nt)] - z[i], 2) - pow(theta[i] - beta[i % (Ng * Nt)] - z[i], 2));
		if (exp(r_t) > R::runif(0, 1)) theta[i] = theta_star;
	}
}
//[[Rcpp::export]]
void update_beta(arma::mat& beta,
	arma::cube theta, arma::rowvec tau2, arma::cube z, // parameters
	arma::mat Sig_b_i, arma::vec eta,                  // hyperparameters
	int Ng, int Nt, int Ns                             // metaparameters
) {
	mat tmz = sum(theta - z, 2);
	for (int i = 0; i < Ng * Nt; i++) {
		double Sig_b_star = 1 / (Ns / tau2[i % Ng] + Sig_b_i.diag()[i]);
		double xi_star = Sig_b_star * (tmz[i] / tau2[i % Ng] + Sig_b_i.diag()[i] * eta[i]);
		beta[i] = R::rnorm(xi_star, sqrt(Sig_b_star));
	}
}
//[[Rcpp::export]]
void update_tau2(arma::rowvec& tau2,
	arma::cube theta, arma::mat beta, arma::cube z, // parameters
	double at, double bt,                           // hyperparameters
	int Ng, int Nt, int Ns                          // metaparameters
) {
	double at_star = (Ns * Nt) / 2 + at;
	for (int g = 0; g < Ng; g++) {
		double tbz = 0;
		for (int i = 0; i < Nt * Ns; i++) tbz += pow(theta.row(g)[i] - beta.row(g)[i % Nt] - z.row(g)[i], 2);
		double bt_star = tbz / 2 + bt;
		tau2[g] = 1 / R::rgamma(at_star, 1 / bt_star);
	}
}
//[[Rcpp::export]]
void update_Gt_i(arma::cube& Gt_i,
	arma::cube z, arma::rowvec rho, arma::mat G, // parameters
	double nu, double bt,                        // hyperparameters
	int Ng, int Nt, int Ns, int I,               // metaparameters
	arma::field<arma::uvec> neigh, arma::vec num // adjacency
) {
	cube Psi_star(Ng, Ng, Nt, fill::zeros);
  cube A(Ng, Ng, 4);
  double nu_Gt_star = Ns - I + nu;
  mat r  = rho.t();
  mat sr = sqrt(1 - pow(r, 2));
  for (int i = 0; i < Ns; i++) {
  	mat sum_zjg = sum(z.slices(neigh[i]), 2);
  	mat zpz     = z.slice(i) - sum_zjg / num[i];
  	vec z0      = z.subcube(0, 0, i, Ng - 1, 0, i);
  	Psi_star.slice(0) += num[i] * (zpz.col(0) * z0.t());
  	for (int t = 1; t < Nt; t++) {
  		vec zt0 = z.subcube(0, t    , i, Ng - 1, t    , i);
  		vec zt1 = z.subcube(0, t - 1, i, Ng - 1, t - 1, i);
  		A.slice(0) =  (1 / sr % zpz.col(t    )) * (1 / sr % zt0).t();
  		A.slice(1) =  (r / sr % zpz.col(t - 1)) * (r / sr % zt1).t();
  		A.slice(2) = -(r / sr % zpz.col(t - 1)) * (1 / sr % zt0).t();
  		A.slice(3) = -(1 / sr % zpz.col(t    )) * (r / sr % zt1).t();
  		Psi_star.slice(t) += num[i] * sum(A, 2);
  	}
  }
  for (int t = 0; t < Nt; t++) {
  	Psi_star.slice(t) += G;
  	String slicenum = "Gt_i";
  	slicenum += t;
  	Gt_i.slice(t) = inv(riwish(nu_Gt_star, Psi_star.slice(t)));
  }
}
//[[Rcpp::export]]
void update_G(arma::mat& G,
	arma::cube Gt_i,                        // parameters
	double nu, double nu_0, arma::mat G0_i, // hyperparameters
	int Ng, int Nt                          // metaparameters
) {
	double nu_G_star = nu * Nt + nu_0;
	mat V_star = inv((mat)sum(Gt_i, 2) + G0_i);
	G = rwish(nu_G_star, V_star);
} 
//[[Rcpp::export]]
void update_z(arma::cube& z,
	arma::cube Gt_i, arma::rowvec rho, arma::rowvec tau2, // parameters
	arma::cube theta, arma::mat beta,                     //
	int Ng, int Nt, int Ns,                               // metaparameters
	arma::field<arma::uvec> neigh, arma::vec num          // adjacency
) {
	mat Sei = Sig_eta_i(Gt_i, rho, Nt);
	mat Taukt_mat(Ng * Nt, Ng * Nt, fill::zeros);
	Taukt_mat.diag() = repmat(1 / tau2, 1, Nt);
	for (int i = 0; i < Ns; i++) {
		mat Sig_z_star = Taukt_mat + num[i] * Sei;
		Sig_z_star     = inv(Sig_z_star);
		vec sum_zj     = mat2vec((mat)sum(z.slices(neigh[i]), 2));
		vec phi_star   = Sig_z_star * (Taukt_mat * mat2vec(theta.slice(i) - beta) + Sei * sum_zj);
		z.slice(i) = vec2mat(rmvnorm(1, phi_star, Sig_z_star), Ng, Nt);
	}
	for (int g = 0; g < Ng; g++) {
		for (int t = 0; t < Nt; t++) {
			z.tube(g, t) -= accu(z.tube(g, t)) / (double)(Ns);
		}
	}
}
//[[Rcpp::export]]
void update_rho(arma::rowvec& rho,
	arma::cube Gt_i, arma::cube z,                   // parameters
	double a_rho, double b_rho, arma::rowvec delta,  // hyperparameters
	int Ng, int Nt, int Ns, int I,                   // metaparameters
	arma::field<arma::uvec> neigh, arma::vec num     // adjacency
) {
	for (int k = 0; k < Ng; k++) {
	  rowvec rho_star = rho;
	  rho_star[k]     = expit(R::rnorm(logit(rho[k]), delta[k]));
		mat Sig_eta_ip  = Sig_eta_i(Gt_i, rho_star, Nt) - Sig_eta_i(Gt_i, rho, Nt);
		double ra = log(1 - pow(rho[k], 2)) - log(1 - pow(rho_star[k], 2));
		double rb = 0;
		for (int i = 0; i < Ns; i++) {
			vec zi    = mat2vec(z.slice(i));
			vec szj   = mat2vec(sum(z.slices(neigh[i]), 2) / num[i]);
			vec zmikt = zi - szj;
			rb += (num[i] * (zi.t() * Sig_eta_ip * zmikt))[0];
		}
		rb /= 2;
		double rc  = a_rho * log(rho_star[k] / rho[k]) + b_rho * log((1 - rho_star[k]) / (1 - rho[k]));
		double r_r = (Ns - I) * (Nt - 1) / 2 * ra - rb + rc;
		if (exp(r_r) > R::runif(0, 1)) rho[k] = rho_star[k];
	}
}
