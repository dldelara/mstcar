// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// progress
void progress(int s, int n_iter, int n, int l);
RcppExport SEXP _mstcar_progress(SEXP sSEXP, SEXP n_iterSEXP, SEXP nSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    progress(s, n_iter, n, l);
    return R_NilValue;
END_RCPP
}
// cube2vec
arma::vec cube2vec(arma::cube c);
RcppExport SEXP _mstcar_cube2vec(SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(cube2vec(c));
    return rcpp_result_gen;
END_RCPP
}
// mat2vec
arma::vec mat2vec(arma::mat m);
RcppExport SEXP _mstcar_mat2vec(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(mat2vec(m));
    return rcpp_result_gen;
END_RCPP
}
// vec2mat
arma::mat vec2mat(arma::rowvec v, int Ng, int Nt);
RcppExport SEXP _mstcar_vec2mat(SEXP vSEXP, SEXP NgSEXP, SEXP NtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type Ng(NgSEXP);
    Rcpp::traits::input_parameter< int >::type Nt(NtSEXP);
    rcpp_result_gen = Rcpp::wrap(vec2mat(v, Ng, Nt));
    return rcpp_result_gen;
END_RCPP
}
// Sig_eta_i
arma::mat Sig_eta_i(arma::cube Gt_i, arma::rowvec rho, int Nt);
RcppExport SEXP _mstcar_Sig_eta_i(SEXP Gt_iSEXP, SEXP rhoSEXP, SEXP NtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Gt_i(Gt_iSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type Nt(NtSEXP);
    rcpp_result_gen = Rcpp::wrap(Sig_eta_i(Gt_i, rho, Nt));
    return rcpp_result_gen;
END_RCPP
}
// acpt_cube
arma::cube acpt_cube(arma::field<arma::cube> f, int Ng, int Nt, int Ns);
RcppExport SEXP _mstcar_acpt_cube(SEXP fSEXP, SEXP NgSEXP, SEXP NtSEXP, SEXP NsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::cube> >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type Ng(NgSEXP);
    Rcpp::traits::input_parameter< int >::type Nt(NtSEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    rcpp_result_gen = Rcpp::wrap(acpt_cube(f, Ng, Nt, Ns));
    return rcpp_result_gen;
END_RCPP
}
// acpt_vec
arma::rowvec acpt_vec(arma::field<arma::mat> f, int Ng);
RcppExport SEXP _mstcar_acpt_vec(SEXP fSEXP, SEXP NgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type Ng(NgSEXP);
    rcpp_result_gen = Rcpp::wrap(acpt_vec(f, Ng));
    return rcpp_result_gen;
END_RCPP
}
// get_outname
String get_outname(String name, String dir, String param, int iter);
RcppExport SEXP _mstcar_get_outname(SEXP nameSEXP, SEXP dirSEXP, SEXP paramSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type name(nameSEXP);
    Rcpp::traits::input_parameter< String >::type dir(dirSEXP);
    Rcpp::traits::input_parameter< String >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(get_outname(name, dir, param, iter));
    return rcpp_result_gen;
END_RCPP
}
// sym_test
arma::mat sym_test(arma::mat A);
RcppExport SEXP _mstcar_sym_test(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(sym_test(A));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_sampler
void gibbs_sampler(List mod, int n_iter, int n_loop, int l);
RcppExport SEXP _mstcar_gibbs_sampler(SEXP modSEXP, SEXP n_iterSEXP, SEXP n_loopSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mod(modSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< int >::type n_loop(n_loopSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    gibbs_sampler(mod, n_iter, n_loop, l);
    return R_NilValue;
END_RCPP
}
// output_cube
arma::field<arma::cube> output_cube(List mod, String param, int burn, int thin, arma::vec file_suff);
RcppExport SEXP _mstcar_output_cube(SEXP modSEXP, SEXP paramSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP file_suffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mod(modSEXP);
    Rcpp::traits::input_parameter< String >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type file_suff(file_suffSEXP);
    rcpp_result_gen = Rcpp::wrap(output_cube(mod, param, burn, thin, file_suff));
    return rcpp_result_gen;
END_RCPP
}
// output_mat
arma::field<arma::mat> output_mat(List mod, String param, int burn, int thin, arma::vec file_suff);
RcppExport SEXP _mstcar_output_mat(SEXP modSEXP, SEXP paramSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP file_suffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mod(modSEXP);
    Rcpp::traits::input_parameter< String >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type file_suff(file_suffSEXP);
    rcpp_result_gen = Rcpp::wrap(output_mat(mod, param, burn, thin, file_suff));
    return rcpp_result_gen;
END_RCPP
}
// get_output
List get_output(List mod, int burn, int thin, StringVector params, arma::vec file_suff);
RcppExport SEXP _mstcar_get_output(SEXP modSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP paramsSEXP, SEXP file_suffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mod(modSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< StringVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type file_suff(file_suffSEXP);
    rcpp_result_gen = Rcpp::wrap(get_output(mod, burn, thin, params, file_suff));
    return rcpp_result_gen;
END_RCPP
}
// acceptance_ratio_cube
arma::cube acceptance_ratio_cube(List mod, arma::vec file_suff, int burn);
RcppExport SEXP _mstcar_acceptance_ratio_cube(SEXP modSEXP, SEXP file_suffSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mod(modSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type file_suff(file_suffSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(acceptance_ratio_cube(mod, file_suff, burn));
    return rcpp_result_gen;
END_RCPP
}
// acceptance_ratio_mat
arma::rowvec acceptance_ratio_mat(List mod, arma::vec file_suff, int burn);
RcppExport SEXP _mstcar_acceptance_ratio_mat(SEXP modSEXP, SEXP file_suffSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mod(modSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type file_suff(file_suffSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(acceptance_ratio_mat(mod, file_suff, burn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mstcar_progress", (DL_FUNC) &_mstcar_progress, 4},
    {"_mstcar_cube2vec", (DL_FUNC) &_mstcar_cube2vec, 1},
    {"_mstcar_mat2vec", (DL_FUNC) &_mstcar_mat2vec, 1},
    {"_mstcar_vec2mat", (DL_FUNC) &_mstcar_vec2mat, 3},
    {"_mstcar_Sig_eta_i", (DL_FUNC) &_mstcar_Sig_eta_i, 3},
    {"_mstcar_acpt_cube", (DL_FUNC) &_mstcar_acpt_cube, 4},
    {"_mstcar_acpt_vec", (DL_FUNC) &_mstcar_acpt_vec, 2},
    {"_mstcar_get_outname", (DL_FUNC) &_mstcar_get_outname, 4},
    {"_mstcar_sym_test", (DL_FUNC) &_mstcar_sym_test, 1},
    {"_mstcar_gibbs_sampler", (DL_FUNC) &_mstcar_gibbs_sampler, 4},
    {"_mstcar_output_cube", (DL_FUNC) &_mstcar_output_cube, 5},
    {"_mstcar_output_mat", (DL_FUNC) &_mstcar_output_mat, 5},
    {"_mstcar_get_output", (DL_FUNC) &_mstcar_get_output, 5},
    {"_mstcar_acceptance_ratio_cube", (DL_FUNC) &_mstcar_acceptance_ratio_cube, 3},
    {"_mstcar_acceptance_ratio_mat", (DL_FUNC) &_mstcar_acceptance_ratio_mat, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mstcar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
