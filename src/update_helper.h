#ifndef UPDATEHELPER_H
#define UPDATEHELPER_H

#include <RcppArmadillo.h>
#include <RcppDist.h>
arma::vec cube2vec(arma::cube c);
arma::vec mat2vec(arma::mat m);
arma::mat vec2mat(arma::rowvec v, int Ng, int Nt);
double logit(double x);
double expit(double p);
arma::mat Sig_eta_i(arma::cube Gt_i, arma::rowvec rho, int Nt);
Rcpp::String get_outname(Rcpp::String name, Rcpp::String dir, Rcpp::String param, int iter);
void sym_test(arma::mat& A);

#endif