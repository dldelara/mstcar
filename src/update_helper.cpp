#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//use templates on these two functions
//[[Rcpp::export]]
arma::vec cube2vec(arma::cube c) {
	vec v(c.n_elem, fill::zeros);
	for (unsigned int i = 0; i < c.n_elem; i++) v[i] = c[i];
	return v;
}
//[[Rcpp::export]]
arma::vec mat2vec(arma::mat m) {
	vec v(m.n_elem, fill::zeros);
	for (unsigned int i = 0; i < m.n_elem; i++) v[i] = m[i];
	return v;
}
//[[Rcpp::export]]
arma::mat vec2mat(arma::rowvec v, int Ng, int Nt) {
	mat m(Ng, Nt, fill::zeros);
	for (unsigned int i = 0; i < v.n_elem; i++) m[i] = v[i];
	return m;
}
double logit(double x) {
	double lx = log(x / (1 - x));
	return lx;
}
double expit(double p) {
	double ep = exp(p) / (1 + exp(p));
	return ep;
}
//[[Rcpp::export]]
arma::mat Sig_eta_i(arma::cube Gt_i, arma::rowvec rho, int Nt) {
	int Ng = rho.n_elem;
	mat r  = repmat(rho, Ng, 1);
	mat sr = sqrt(1 - pow(r, 2));
	mat Sei(Ng * Nt, Ng * Nt, fill::zeros);
	mat Gti;
	int it1;
  int it2;
  int il1;
  int il2;
	Sei.submat(0, 0, Ng - 1, Ng - 1) = Gt_i.slice(0);
	for (int t = 1; t < Nt; t ++) {
		it1 = t * Ng;
		it2 = (t + 1) * Ng - 1;
		il1 = (t - 1) * Ng;
		il2 = t * Ng - 1;
		Gti = Gt_i.slice(t);
		Sei.submat(il1, il1, il2, il2) += (r  / sr).t() % (r / sr % Gti);
		Sei.submat(it1, it1, it2, it2)  = (1  / sr).t() % (1 / sr % Gti);
		Sei.submat(il1, it1, il2, it2)  = (-r / sr).t() % (1 / sr % Gti);
		Sei.submat(it1, il1, it2, il2)  = (-r / sr) % (1 / sr % Gti).t();
	}
	return Sei;
}
//[[Rcpp::export]]
arma::mat sym_test(arma::mat A) {
	int p = A.n_rows;
	for (int i = 0; i < (p - 1); i++) {
		for (int j = i + 1; j < p; j++) {
			if (abs(A(i, j) - A(j, i)) / A(i, j) > 1e-5) {
				Rcout << A << "\n" << abs(A(i, j) - A(j, i)) / A(i, j) * 100 << "\n";
				stop("Sad :(");
			} 
		}
	}
	return (A + A.t()) / 2;
}
