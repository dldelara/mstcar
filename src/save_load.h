#ifndef SAVELOAD_H
#define SAVELOAD_H

#include <RcppArmadillo.h>
Rcpp::String get_outname(Rcpp::String name, Rcpp::String dir, Rcpp::String param, int iter);
arma::field<arma::cube> output_cube(Rcpp::List mod, Rcpp::String param, int burn, int thin, arma::vec file_suff);
arma::field<arma::mat> output_mat(Rcpp::List mod, Rcpp::String param, int burn, int thin, arma::vec file_suff);
Rcpp::List get_output(Rcpp::List mod, int burn, int thin, Rcpp::StringVector params, arma::vec file_suff);

#endif