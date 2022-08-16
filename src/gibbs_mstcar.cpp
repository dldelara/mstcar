#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void progress(int s, int n_iter, int n, int l) {
  if (s == 0) {
    String progress(" Progress: |..................................................|");
    if (l > 0) Rcout << "Batch: " << n << "/" << l << progress.get_cstring() <<  "\r";
    else       Rcout << progress.get_cstring() << "\r";
  }
  bool prog = false;
  if ((s + 1) % (n_iter / 50) == 0) prog = true;
  if (prog) {
    String progress(" Progress: |..................................................|");
    for (int i = 0; i < round((s + 1) * 50.0 / n_iter); i++) progress.replace_first(".", "*");
    if (l > 0) Rcout << "Batch: " << n << "/" << l << progress.get_cstring() <<  "\r";
    else       Rcout << progress.get_cstring() << "\r";
  }
  if ((s == n_iter - 1) & (n == l)) Rcout << "\n";
}
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
String get_outname(String name, String dir, String param, int iter) {
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
//[[Rcpp::export]]
void gibbs_sampler(List mod, int n_iter, int n_loop = 0, int l = 0) {
	// Import data, parameters, and priors
	List data     = mod["data"];
	List nb       = mod["nb"];
	List params   = mod["params"];
	List priors   = mod["priors"];
	List inits    = mod["inits"];
	cube Y			  = data["Y"];
	cube n        = data["n"];
	vec  dNd      = params["dNd"];
	String name   = params["name"];
	String dir    = params["dir"];
	int  Ng       = dNd[0];
	int  Nt       = dNd[1];
	int  Ns       = dNd[2];
	int  Nd       = Ng * Nt * Ns;
	int  I        = params["I"];
	bool rho_up   = params["rho_up"];
	String method = params["method"];
	int its       = params["its"];
	vec  num      = nb["num"];
	field<uvec> neigh = nb["neigh"];
	vec  eta      = priors["eta"];
	mat  Sig_b_i  = priors["Sig_b_i"];
	mat  G0_i     = priors["G0_i"];
	cube gamma    = priors["gamma"];
	double at     = priors["at"];
	double bt     = priors["bt"];
	double nu     = priors["nu"];
	double nu_0   = priors["nu_0"];
	double a_rho  = 0;
	double b_rho  = 0;
	rowvec delta  = rowvec(Ng, fill::zeros);
	if (rho_up) {
		a_rho = priors["a_rho"];
		b_rho = priors["b_rho"];
		delta = as<rowvec>(priors["delta"]);
	}
  // Import initial values
	cube theta  = inits["theta"];
	cube z      = inits["z"];
	mat  beta   = inits["beta"];
	mat  G      = inits["G"];
	cube Gt_i   = inits["Gt"];
	rowvec rho  = inits["rho"];
	rowvec tau2 = inits["tau2"];
	for (int t = 0; t < Nt; t++) Gt_i.slice(t) = inv(Gt_i.slice(t));
	// Set up output
	field<cube> new_theta(n_iter);
	field<cube> new_Gt   (n_iter);
	field<cube> new_z    (n_iter);
	field<mat>  new_beta (n_iter);
	field<mat>  new_G    (n_iter);
	field<mat>  new_tau2 (n_iter);
	field<mat>  new_rho  (n_iter);
	cube Gt(Ng, Ng, Nt, fill::zeros);
	// Establish sampler parameters
	double theta_star; // theta
	double r_t;
	double Sig_b_star; // beta (diagonal)
	double xi_star;
	mat tmz;
	double at_star = (Ns * Nt) / 2 + at; // tau2
	double bt_star;
	double tbz;
	double nu_G_star; // G
	mat  V_star(Ng, Ng, fill::zeros);
	mat  Sei(Ng * Nt, Ng * Nt, fill::zeros); // z
	mat  Taukt_mat(Ng * Nt, Ng * Nt, fill::zeros);
	mat  Sig_z_star;
	vec  sum_zj;
	vec  phi_star;
	cube Psi_star(Ng, Ng, Nt, fill::zeros); // Gt_i
	double nu_Gt_star;
	mat  sum_zjg;
  mat  r;
	mat  sr;
  mat  zpz;
  vec  zpz0;
  vec  z0;
  vec  zt0;
  vec  zt1;
  cube A(Ng, Ng, 4, fill::zeros);
	cube gam_acpt; // Adaptive variance
	field<cube> theta_a;
	rowvec rho_star; // rho
	double ra;
	double rb;
	double rc;
	double r_r;
	mat Sig_eta_ip;
	vec zmikt;
	vec zi;
	vec szj;
	rowvec rho_acpt;
	field<mat> rho_a;
	for (int s = 0; s < n_iter; s++) {
		// Update theta
		for (int i = 0; i < Nd; i++) {
			theta_star = R::rnorm(theta[i], gamma[i]);
			r_t = Y[i] * (theta_star - theta[i]);
			if (method == "binom") r_t += - n[i] * (log(1 + exp(theta_star)) - log(1 + exp(theta[i])));
			if (method == "pois" ) r_t += - n[i] * (exp(theta_star) - exp(theta[i]));
			r_t += - 1 / (2 * tau2[i % Ng]) * (pow(theta_star - beta[i % (Ng * Nt)] - z[i], 2) - pow(theta[i] - beta[i % (Ng * Nt)] - z[i], 2));
			if (exp(r_t) > R::runif(0, 1)) theta[i] = theta_star;
		}
		// Update beta
		// Add contingency for non-diagonal prior
		if (Sig_b_i.is_diagmat()) {
			tmz = sum(theta - z, 2);
			for (int i = 0; i < Ng * Nt; i++) {
				Sig_b_star = 1 / (Ns / tau2[i % Ng] + Sig_b_i.diag()[i]);
				xi_star = Sig_b_star * (tmz[i] / tau2[i % Ng] + Sig_b_i.diag()[i] * eta[i]);
				beta[i] = R::rnorm(xi_star, sqrt(Sig_b_star));
			}
		}
		
		// Update tau2
		for (int g = 0; g < Ng; g++) {
			tbz = 0;
			for (int i = 0; i < Nt * Ns; i++) tbz += pow(theta.row(g)[i] - beta.row(g)[i % Nt] - z.row(g)[i], 2);
			bt_star = tbz / 2 + bt;
			tau2[g] = 1 / R::rgamma(at_star, 1 / bt_star);
		}
		// Update Gt_i
		Psi_star.zeros();
	  nu_Gt_star = Ns - I + nu;
	  r  = rho.t();
	  sr = sqrt(1 - pow(r, 2));
	  for (int i = 0; i < Ns; i++) {
	  	sum_zjg = sum(z.slices(neigh[i]), 2);
	  	zpz     = z.slice(i) - sum_zjg / num[i];
	  	z0      = z.subcube(0, 0, i, Ng - 1, 0, i);
	  	Psi_star.slice(0) += num[i] * (zpz.col(0) * z0.t());
	  	for (int t = 1; t < Nt; t++) {
	  		zt0 = z.subcube(0, t    , i, Ng - 1, t    , i);
	  		zt1 = z.subcube(0, t - 1, i, Ng - 1, t - 1, i);
	  		A.slice(0) =  (1 / sr % zpz.col(t    )) * (1 / sr % zt0).t();
	  		A.slice(1) =  (r / sr % zpz.col(t - 1)) * (r / sr % zt1).t();
	  		A.slice(2) = -(r / sr % zpz.col(t - 1)) * (1 / sr % zt0).t();
	  		A.slice(3) = -(1 / sr % zpz.col(t    )) * (r / sr % zt1).t();
	  		Psi_star.slice(t) += num[i] * sum(A, 2);
	  	}
	  }
	  for (int t = 0; t < Nt; t++) {
	  	Psi_star.slice(t) += G;
	  	if (!(Psi_star.slice(t).is_symmetric())) Psi_star.slice(t) = sym_test(Psi_star.slice(t));
	  	Gt_i.slice(t) = inv(riwish(nu_Gt_star, Psi_star.slice(t)));
	  }
		// Update G
		nu_G_star = nu * Nt + nu_0;
		V_star = sum(Gt_i, 2);
		V_star = inv(V_star + G0_i);
		G      = rwish(nu_G_star, V_star);
		// Update z
		Sei = Sig_eta_i(Gt_i, rho, Nt);
		Taukt_mat.diag() = repmat(1 / tau2, 1, Nt);
		for (int i = 0; i < Ns; i++) {
			Sig_z_star = Taukt_mat + num[i] * Sei;
			if (!(Sig_z_star.is_symmetric())) Sig_z_star = sym_test(Sig_z_star);
			Sig_z_star = inv(Sig_z_star);
			sum_zj     = cube2vec(sum(z.slices(neigh[i]), 2));
			phi_star   = Sig_z_star * (Taukt_mat * mat2vec(theta.slice(i) - beta) + Sei * sum_zj);
			z.slice(i) = vec2mat(rmvnorm(1, phi_star, Sig_z_star), Ng, Nt);
		}
		for (int g = 0; g < Ng; g++) {
			for (int t = 0; t < Nt; t++) {
				z.tube(g, t) -= accu(z.tube(g, t)) / (double)(Ns);
			}
		}
		// Update rho
		if (rho_up) {
			for (int k = 0; k < Ng; k++) {
			  rho_star = rho;
			  rho_star[k] = expit(R::rnorm(logit(rho[k]), delta[k]));
				ra = log(1 - pow(rho[k], 2)) - log(1 - pow(rho_star[k], 2));
				Sig_eta_ip = Sig_eta_i(Gt_i, rho_star, Nt) - Sig_eta_i(Gt_i, rho, Nt);
				rb = 0;
				for (int i = 0; i < Ns; i++) {
					zi    = mat2vec(z.slice(i));
					szj   = mat2vec(sum(z.slices(neigh[i]), 2) / num[i]);
					zmikt = zi - szj;
					rb   += (num[i] * (zi.t() * Sig_eta_ip * zmikt))[0];
				}
				rb /= 2;
				rc  = a_rho * log(rho_star[k] / rho[k]) + b_rho * log((1 - rho_star[k]) / (1 - rho[k]));
				r_r = (Ns - I) * (Nt - 1) / 2 * ra - rb + rc;
				if (exp(r_r) > R::runif(0, 1)) rho[k] = rho_star[k];
			}
		}
		// Adaptive variance updates
		if (((s + 1.0) / 100 == floor((s + 1.0) / 100)) | (s == n_iter - 1)) {
			theta_a  = new_theta.rows(s - 99, s);
			gam_acpt = acpt_cube(theta_a, Ng, Nt, Ns);
			for (int i = 0; i < Nd; i++) {
				if (gam_acpt[i] > 0.75) gam_acpt[i] = 0.75;
				if (gam_acpt[i] < 0.2) gam_acpt[i] = 0.2;
				gamma[i] *= gam_acpt[i] / 0.43;
			}
			if (rho_up) {
				rho_a = new_rho.rows(s - 99, s);
				rho_acpt = acpt_vec(rho_a, Ng);
				for (int k = 0; k < Ng; k++) {
					if (rho_acpt[k] > 0.75) rho_acpt[k] = 0.75;
					if (rho_acpt[k] < 0.2) rho_acpt[k] = 0.2;
					delta[k] *= rho_acpt[k] / 0.43;
				}
			}
		}

		for (int t = 0; t < Nt; t++) Gt.slice(t) = inv(Gt_i.slice(t));
		// Store samples
		new_theta[s] = theta;
		new_beta [s] = beta;
		new_z    [s] = z;
		new_tau2 [s] = tau2;
		new_Gt   [s] = Gt;
		new_G    [s] = G;
		if (rho_up) new_rho[s] = rho;
		progress(s, n_iter, n_loop, l);
	}
	// Save new inits for next samples
	mod["inits"] = List::create(
		_["theta"] = theta,
		_["beta"]  = beta,
		_["tau2"]  = tau2,
		_["Gt"]    = Gt,
		_["G"]     = G,
		_["z"]     = z,
		_["rho"]   = rho
	);
	// Update parameters and priors
	params["its"] = its + n_iter;
	mod["params"] = params;
	if (rho_up) priors["delta"] = delta;
	priors["gamma"] = gamma;
	mod["priors"]   = priors;
	new_theta.save(get_outname(name, dir, "theta", its + n_iter).get_cstring());
	new_beta .save(get_outname(name, dir, "beta" , its + n_iter).get_cstring());
	new_tau2 .save(get_outname(name, dir, "tau2" , its + n_iter).get_cstring());
	new_G    .save(get_outname(name, dir, "G"    , its + n_iter).get_cstring());
	new_Gt   .save(get_outname(name, dir, "Gt"   , its + n_iter).get_cstring());
	new_z    .save(get_outname(name, dir, "z"    , its + n_iter).get_cstring());
	if (rho_up) new_rho.save(get_outname(name, dir, "rho", its + n_iter).get_cstring());
	field<mat> output_test;
	output_test.load(get_outname(name, dir, "tau2", its + n_iter).get_cstring());
	for (int i = 0; i < 2; i++) {
		if (size(output_test)[i] == 0) {
			stop("Output saved improperly. Please try re-creating model and running again.");
		}
	}
}
//[[Rcpp::export]]
arma::field<arma::cube> output_cube(List mod, String param, int burn, int thin, arma::vec file_suff) {
	List params   = mod["params"];
	int its       = params["its"];
	String name   = params["name"];
	String dir    = params["dir"];
	String method = params["method"];
	String file;
	field<cube> output_full;
	field<cube> output_thin((its - burn) / thin);
	int j = 0;
	for (unsigned int it = 0; it < file_suff.n_elem; it++) {
		file = get_outname(name, dir, param.get_cstring(), file_suff[it]);
		Rcout << "Pulling files from: " << file.get_cstring() << "\n"; 
		output_full.load(file.get_cstring());
		for (unsigned int i = thin - 1; i < output_full.n_elem; i += thin) {
			if ((param == "theta") | (param == "z")) {
				if (method == "binom") {
					output_thin[j] = exp(output_full[i]) / (1 + exp(output_full[i]));
				}
				if (method == "pois") {
					output_thin[j] = exp(output_full[i]);
				}
			}
			else output_thin[j] = output_full[i];
			j++;
		}
	}
	return output_thin;
}
//[[Rcpp::export]]
arma::field<arma::mat> output_mat(List mod, String param, int burn, int thin, arma::vec file_suff) {
	List params   = mod["params"];
	int its       = params["its"];
	String name   = params["name"];
	String dir    = params["dir"];
	String method = params["method"];
	String file;
	field<mat> output_full;
	field<mat> output_thin((its - burn) / thin);
	int j = 0;
	for (unsigned int it = 0; it < file_suff.n_elem; it++) {
		file = get_outname(name, dir, param.get_cstring(), file_suff[it]);
		Rcout << "Pulling files from: " << file.get_cstring() << "\n";
		output_full.load(file.get_cstring());
		for (unsigned int i = thin - 1; i < output_full.n_elem; i += thin) {
			if (param == "beta") {
				if (method == "binom") {
					output_thin[j] = exp(output_full[i]) / (1 + exp(output_full[i]));
				}
				if (method == "pois") {
					output_thin[j] = exp(output_full[i]);
				}
			}
			else output_thin[j] = output_full[i];
			j++;
		}
	}
	return output_thin;
}
//[[Rcpp::export(".load_samples")]]
List load_samples(List mod, StringVector params, int burn, int thin, arma::vec file_suff) {
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
//[[Rcpp::export]]
arma::cube acceptance_ratio_cube(List mod, arma::vec file_suff, int burn) {
	List params = mod["params"];
	String name = params["name"];
	String dir  = params["dir"];
	vec dNd     = params["dNd"];
	int Ng      = dNd[0];
	int Nt      = dNd[1];
	int Ns      = dNd[2];
	List sample = load_samples(mod, "theta", burn, 1, file_suff);
	field<cube> theta = sample["theta"];
	return acpt_cube(theta, Ng, Nt, Ns);
}
//[[Rcpp::export]]
arma::rowvec acceptance_ratio_mat(List mod, arma::vec file_suff, int burn) {
	List params = mod["params"];
	String name = params["name"];
	String dir  = params["dir"];
	vec dNd     = params["dNd"];
	int Ng      = dNd[0];
	List sample = load_samples(mod, "rho", burn, 1, file_suff);
	field<mat> rho = sample["rho"];
	return acpt_vec(rho, Ng);
}
