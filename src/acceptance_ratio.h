#ifndef ACCEPTANCERATIO_H
#define ACCEPTANCERATIO_H

#include <RcppArmadillo.h>
#include <RcppDist.h>
arma::cube acpt_cube(arma::field<arma::cube> f, int Ng, int Nt, int Ns);
arma::rowvec acpt_vec(arma::field<arma::mat> f, int Ng);
arma::cube acceptance_ratio_cube(Rcpp::List mod, arma::vec file_suff, int burn);
arma::rowvec acceptance_ratio_mat(Rcpp::List mod, arma::vec file_suff, int burn);
void adaptive_variance_theta(arma::cube& gamma, arma::field<arma::cube> new_theta, int s, int Ng, int Nt, int Ns);
void adaptive_variance_rho(arma::rowvec& delta, arma::field<arma::mat> new_rho, int s, int Ng);

#endif