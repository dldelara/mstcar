#ifndef PARAMUPDATES_H
#define PARAMUPDATES_H

#include <RcppArmadillo.h>
#include <RcppDist.h>
void update_theta(arma::cube& theta,
	arma::cube Y, arma::cube n,
	arma::rowvec tau2, arma::mat beta, arma::cube z,
	arma::cube gamma,
	int Ng, int Nt, int Ns, Rcpp::String method);
void update_beta(arma::mat& beta,
	arma::cube theta, arma::rowvec tau2, arma::cube z,
	arma::mat Sig_b_i, arma::vec eta,
	int Ng, int Nt, int Ns);
void update_tau2(arma::rowvec& tau2,
	arma::cube theta, arma::mat beta, arma::cube z,
	double at, double bt,
	int Ng, int Nt, int Ns);
void update_Gt_i(arma::cube& Gt_i,
	arma::cube z, arma::rowvec rho, arma::mat G,
	double nu, double bt,
	int Ng, int Nt, int Ns, int I,
	arma::field<arma::uvec> neigh, arma::vec num);
void update_G(arma::mat& G,
	arma::cube Gt_i,
	double nu, double nu_0, arma::mat G0_i,
	int Ng, int Nt);
void update_z(arma::cube& z,
	arma::cube Gt_i, arma::rowvec rho, arma::rowvec tau2,
	arma::cube theta, arma::mat beta,
	int Ng, int Nt, int Ns,
	arma::field<arma::uvec> neigh, arma::vec num);
void update_rho(arma::rowvec& rho,
	arma::cube Gt_i, arma::cube z,
	double a_rho, double b_rho, arma::rowvec delta,
	int Ng, int Nt, int Ns, int I,
	arma::field<arma::uvec> neigh, arma::vec num);

#endif