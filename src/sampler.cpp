#include <RcppArmadillo.h>
#include "param_updates.h"
#include "acceptance_ratio.h"
#include "progress.h"
#include "save_load.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void gibbs_sampler(List& mod, int n_iter, int n_loop = 0, int l = 0) {
	// Import data, parameters, and priors
	List data     = mod["data"];
	List nb       = mod["nb"];
	List params   = mod["params"];
	List priors   = mod["priors"];
	List inits    = mod["inits"];
	cube Y			  = data["Y"];
	cube n        = data["n"];
	String name   = params["name"];
	String dir    = params["dir"];
	bool rho_up   = params["rho_up"];
	int  its      = params["its"];
	int  I        = params["I"];
	String method = params["method"];
	vec dNd       = params["dNd"];
	int Ng        = dNd[0];
	int Nt        = dNd[1];
	int Ns        = dNd[2];
	field<uvec> neigh = nb["neigh"];
	vec  num      = nb["num"];
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
	// Gibbs sampler
	progress(0, n_iter, n_loop, l);
	for (int s = 0; s < n_iter; s++) {
		update_theta(theta,  //// Update theta
			Y, n,              // data
		  tau2, beta, z,     // parameters
		  gamma,             // hyperparameters
		  Ng, Nt, Ns, method // metaparameters
		);
		update_beta(beta, //// Update beta
			theta, tau2, z, // parameters
			Sig_b_i, eta,   // hyperparameters
			Ng, Nt, Ns      // metaparameters
		);
		update_tau2(tau2, //// Update tau2
			theta, beta, z, // parameters
			at, bt,         // hyperparameters
			Ng, Nt, Ns      // metaparameters
		);
		update_Gt_i(Gt_i, //// Update Gt_i
			z, rho, G,      // parameters
			nu, bt,         // hyperparameters
			Ng, Nt, Ns, I,  // metaparameters
			neigh, num      // adjacency
		);
		update_G(G,       //// Update G
			Gt_i,           // parameters
			nu, nu_0, G0_i, // hyperparameters
			Ng, Nt          // metaparameters
		);
		update_z(z,        //// Update z
    	Gt_i, rho, tau2, // parameters
    	theta, beta,     //
    	Ng, Nt, Ns,      // metaparameters
    	neigh, num       // adjacency
    );
		if (rho_up) 
			update_rho(rho,        //// Update rho
		 		Gt_i, z,             // parameters
				a_rho, b_rho, delta, // hyperparameters
				Ng, Nt, Ns, I,       // metaparameters
				neigh, num           // adjacency
			);
		// Adaptive variance updates
		if ((s + 1) % 100 == 0) {
			adaptive_variance_theta(gamma, new_theta, s, Ng, Nt, Ns);
			if (rho_up) adaptive_variance_rho(delta, new_rho, s, Ng);
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
	new_Gt   .save(get_outname(name, dir, "Gt"   , its + n_iter).get_cstring());
	new_G    .save(get_outname(name, dir, "G"    , its + n_iter).get_cstring());
	new_z    .save(get_outname(name, dir, "z"    , its + n_iter).get_cstring());
	if (rho_up) new_rho.save(get_outname(name, dir, "rho", its + n_iter).get_cstring());
}

