#' Calculate Acceptance Ratio for Metropolis Parameters
#'
#' @param mod a model object created with the \code{init()} function.
#' @param params a string vector of parameters to load samples for. \code{all} specifies \code{theta} and \code{rho}.
#' @param burn amount of burn-in to be done on samples. Discards first \code{burn} iterations of the samples.
#' @return A list containing the acceptance ratios for each \code{theta} and/or \code{rho} parameter.
#' @examples
#' \dontrun{
#' acceptance_ratio(mod_nc)
#' acceptance_ratio(mod_nc, params = "all")
#' acceptance_ratio(mod = mod_nc, params = "theta", burn = 2e3)
#' }
#' @export
acceptance_ratio = function(mod, params = c("all", "theta", "rho"), burn = 0) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all" %in% params)  params = c("theta", "rho")
  if (!mod$params$rho_up) params = params[which(params != "rho")]
  accept = NULL
  if ("theta" %in% params) accept$theta = acceptance_ratio_cube(mod, getsuff(mod, "theta", burn), burn)
  if ("rho"   %in% params) accept$rho   = acceptance_ratio_mat (mod, getsuff(mod, "rho"  , burn), burn)
  if (any(c(accept) > 0.47) | any(c(accept) < 0.23)) warning("MCMC algorithm is still in its adaption phase. More iterations are required.\n")
  accept
}
#' Calculate reliability of median estimates
#'
#' @param medians a list of medians created with the \code{get_medians()} function.
#' @param params a string vector of parameters to load samples for. \code{all} specifies all parameters.
#' @return A list containing median estimates for each parameter of interest.
#' @examples
#' \dontrun{
#' get_reliability(medians_nc)
#' get_reliability(medians_nc, c("theta", "tau2", "Gt"))
#' }
#' @export
get_reliability = function(medians, params = c("all", names(medians))) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all" %in% params) params = names(medians)
  rel = lapply(medians[params], \(x) apply(x, 2:length(dim(x)), \(.) .[1] > (.[3] - .[2])))
  rel
}
