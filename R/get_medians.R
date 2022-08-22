#' Calculate median estimates from samples
#'
#' @param rec_samples a list of samples created with the \code{load_samples()} function.
#' @param params a string vector of parameters to load samples for. \code{all} specifies all parameters.
#' @param ci the credible interval to be calculated for estimates.
#' @return A list containing median estimates for each parameter of interest.
#' @examples
#' \dontrun{
#' load_samples(output_nc)
#' load_samples(output_nc, params = "all")
#' load_samples(output_nc, c("theta", "tau2", "Gt"))
#' }
#' @export
get_medians = function(rec_samples, params = c("all", names(rec_samples)), ci = 0.95) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all"  %in% params) params = names(rec_samples)
  if ("tau2" %in% params) rec_samples$tau2 = t(rec_samples$tau2)
  if ("rho"  %in% params) rec_samples$rho  = t(rec_samples$rho )
  est = NULL
  for (nam in params) {
    d = 1:(length(dim(rec_samples[[nam]])) - 1)
    est[[nam]] = apply(rec_samples[[nam]], d, quantile, c(0.5, (1 - ci) / 2, 1 - (1 - ci) / 2))
  }
  est
}
