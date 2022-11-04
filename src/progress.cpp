#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void progress(int s, int n_iter, int n, int l) {
  if (s == 0) {
    String progress(" Progress: |..................................................|");
    if (l > 0) Rcout << "Batch: " << n << "/" << l << progress.get_cstring() << "\r";
    else       Rcout << progress.get_cstring() << "\r";
  }
  bool prog = false;
  if ((s + 1) % (n_iter / 50) == 0) prog = true;
  if (prog) {
    String progress(" Progress: |..................................................|");
    for (int i = 0; i < round((s + 1) * 50.0 / n_iter); i++) progress.replace_first(".", "*");
    if (l > 0) Rcout << "Batch: " << n << "/" << l << progress.get_cstring() << "\r";
    else       Rcout << progress.get_cstring() << "\r";
  }
  if ((s == n_iter - 1) & (n == l)) Rcout << "\n";
}