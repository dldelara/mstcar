#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::String get_outname(Rcpp::String name, Rcpp::String dir, Rcpp::String param, int iter) {
	String out_string;
	out_string.push_back(dir);
	out_string.push_back("/");
	out_string.push_back(name);
	out_string.push_back("/");
	out_string.push_back(param);
	out_string.push_back("/");
	out_string.push_back(name);
	out_string.push_back("_");
	out_string.push_back(param);
	out_string.push_back("_");
	out_string.push_back(iter);
	out_string.push_back(".txt");
	return out_string;
}
// use templates on these two functions
//[[Rcpp::export]]
arma::field<arma::cube> output_cube(Rcpp::List mod, Rcpp::String param, int burn, int thin, arma::vec file_suff) {
	List params   = mod["params"];
	int its       = params["its"];
	String name   = params["name"];
	String dir    = params["dir"];
	String method = params["method"];
	field<cube> output_full;
	field<cube> output_thin((its - burn) / thin);
	for (unsigned int it = 0, j = 0; it < file_suff.n_elem; it++) {
		String file = get_outname(name, dir, param.get_cstring(), file_suff[it]);
		output_full.load(file.get_cstring());
		for (unsigned int i = thin - 1; i < output_full.n_elem; i += thin, j++) {
			if ((param == "theta") | (param == "z")) {
				if (method == "binom") output_thin[j] = exp(output_full[i]) / (1 + exp(output_full[i]));
				if (method == "pois")  output_thin[j] = exp(output_full[i]);
			} else output_thin[j] = output_full[i];
		}
	}
	return output_thin;
}
//[[Rcpp::export]]
arma::field<arma::mat> output_mat(Rcpp::List mod, Rcpp::String param, int burn, int thin, arma::vec file_suff) {
	List params   = mod["params"];
	int its       = params["its"];
	String name   = params["name"];
	String dir    = params["dir"];
	String method = params["method"];
	field<mat> output_full;
	field<mat> output_thin((its - burn) / thin);
	for (unsigned int it = 0, j = 0; it < file_suff.n_elem; it++) {
		String file = get_outname(name, dir, param.get_cstring(), file_suff[it]);
		output_full.load(file.get_cstring());
		for (unsigned int i = thin - 1; i < output_full.n_elem; i += thin, j++) {
			if (param == "beta") {
				if (method == "binom") output_thin[j] = exp(output_full[i]) / (1 + exp(output_full[i]));
				if (method == "pois")  output_thin[j] = exp(output_full[i]);
			} else output_thin[j] = output_full[i];
		}
	}
	return output_thin;
}
//[[Rcpp::export]]
Rcpp::List get_output(Rcpp::List mod, int burn, int thin, Rcpp::StringVector params, arma::vec file_suff) {
	List samples;
	for (int p = 0; p < params.length(); p++) {
		String param = params[p];
		if ((param == "theta") | (param == "Gt") | (param == "z")) {
			samples[param.get_cstring()] = output_cube(mod, param.get_cstring(), burn, thin, file_suff);		
		} else {
			samples[param.get_cstring()] = output_mat (mod, param.get_cstring(), burn, thin, file_suff);	
		} 
	}
	return samples;
}
