#' Initialize MSTCAR model for Gibbs sampling
#'
#' @param data List of binomial- or poisson-distributed data. Requires a list containing 2 objects:
#' \code{n} and \code{Y}. \code{n} represents total population and \code{Y}
#' represents event data.
#' @param nb Optional list of adjacency data, in WinBugs format. Requires a list containing 2 objects:
#' \code{neigh} and \code{num}. \code{neigh} is a list of length \code{Ns},
#' the number of regions, describing the neighbors for each region. \code{num}
#' is a vector of length \code{Ns} marking the number of neighbors in each region.
#' \code{neigh} is required if \code{shp} is not defined in \code{data}.
#' @param shp A dataframe object that includes a variable called \code{geometry}. \code{geometry}
#' observations include a list of each coordinate in the shapefile. This is required
#' if \code{nb} is not specified. Every region in \code{shp} must have at least one neighbor.
#' @param priors Optional list of hyperparameters. If unspecified, will use default
#' values. Requires a list containing at least 1 of 8 optional objects: \code{gamma}, \code{eta},
#' \code{Sig_b_i}, \code{at}, \code{bt}, \code{nu}, \code{nu_0}, and \code{G0}. If \code{rho_up = TRUE},
#' will also take \code{a_rho}, \code{b_rho}, and \code{delta}.
#' @param inits Optional list of initial values for Gibbs sampler.
#' If unspecified, will use default values. Requires a list containing at least 1 of 7 objects:
#' \code{theta}, \code{beta}, \code{tau2}, \code{Gt}, \code{G}, \code{z}, and \code{rho}.
#' @param name the name given to the model. This will decide the name of the directory that's created
#' for the model output.
#' @param dir the directory that the model will be saved to. Defaults to \code{getwd()}.
#' @param rho_up decides whether \code{rho} will be updated or treated as a constant.
#' @param method decides which method to use to calculate samples for the \code{theta} Metropolis
#' updates: \code{method = "binom"} will do Binomial updates and \code{method = "pois"} will do
#' Poisson updates.
#' @return A list containing model data, neighbor data, priors, initial output,
#' and auxiliary hyperparameters related to the data.
#' @examples
#' \dontrun{
#' init(data)
#' init(data = ncheart, nb = ncnb)
#' init(ncheart, shp = ncshp, priors = us_priors)
#' }
#' @export
init = function(
    data,
    nb     = NULL,
    shp    = NULL,
    priors = NULL,
    inits  = NULL,
    name,
    dir    = getwd(),
    rho_up = FALSE,
    method = c("binom", "pois")
) {
  method = match.arg(method)
  mod = list(data = data, nb = nb, shp = shp, priors = priors, inits = inits, params = NULL)
  # Shapefile and neighbor data missing
  if (is.null(mod$shp) & is.null(mod$nb$neigh)) {
    stop("Shapefile and neighbor data missing. Please specify either shapefile, neighbor data, or both")
  }
  nbmiss = NULL
  if (!is.null(mod$shp)) {
    if (is.null(mod$nb$neigh)) {
      mod$nb$neigh = spdep::poly2nb(mod$shp)
      nbmiss       = c(nbmiss, "neigh")
    }
    if (is.null(mod$nb$num)) {
      mod$nb$num = sapply(mod$nb$neigh, length)
      nbmiss     = c(nbmiss, "num")
    }
  }
  if (!is.null(nbmiss)) cat("The following objects were created using defaults in 'nb':", paste(nbmiss, collapse = " "), "\n")
  nbcheck(mod)
  mod$params$method = method
  mod$params$rho_up = rho_up
  mod$params$its    = 0
  mod$params$dNd    = dim(mod$data$Y)
  mod$params$I      = length(getislands(mod))
  mod$params$name   = name
  mod$params$dir    = gsub('(/)\\1+', '\\1', dir)
  mod$nb$neigh      = sapply(mod$nb$neigh, \(x) x - 1) # translate neighbor information to zero-index
  Ng = dim(mod$data$Y)[1]
  Nt = dim(mod$data$Y)[2]
  # Prior parameter check
  primiss = NULL
  if (is.null(mod$priors$eta)) {
    mod$priors$eta = rep(-5, Ng * Nt)
    primiss = c(primiss, "eta")
  }
  if (is.null(mod$priors$Sig_b_i)) {
    mod$priors$Sig_b_i = diag(1e-4, Ng * Nt)
    primiss = c(primiss, "Sig_b_i")
  }
  if (is.null(mod$priors$at)) {
    mod$priors$at = 100
    primiss = c(primiss, "at")
  }
  if (is.null(mod$priors$bt)) {
    mod$priors$bt = 1
    primiss = c(primiss, "bt")
  }
  if (is.null(mod$priors$nu)) {
    mod$priors$nu = Ng + 2
    primiss = c(primiss, "nu")
  }
  if (is.null(mod$priors$nu_0)) {
    mod$priors$nu_0 = Ng + 2
    primiss = c(primiss, "nu_0")
  }
  if (is.null(mod$priors$G0_i)) {
    mod$priors$G0_i = diag(7, Ng)
    primiss = c(primiss, "G0_i")
  }
  if (is.null(mod$priors$gamma)) {
    mod$priors$gamma = array(0.025, dim = mod$params$dNd)
    primiss = c(primiss, "gamma")
  }
  if (rho_up) {
    if (is.null(mod$priors$a_rho)) {
      mod$priors$a_rho = 95
      primiss = c(primiss, "a_rho")
    }
    if (is.null(mod$priors$b_rho)) {
      mod$priors$b_rho = 5
      primiss = c(primiss, "b_rho")
    }
    if (is.null(mod$priors$delta)) {
      mod$priors$delta = rep(0.05, Ng)
      primiss = c(primiss, "delta")
    }
  }

  if (!is.null(primiss)) cat("The following objects were created using defaults in 'priors':", paste(primiss, collapse = " "), "\n")
  pricheck(mod)
  # Initial value check
  initmiss = NULL
  if (is.null(mod$inits$theta)) {
    if (method == "binom") mod$inits$theta = logit(ifelse((data$Y == 0) | (data$n == 0) | (data$Y >= data$n), sum(data$Y) / sum(data$n), data$Y / data$n))
    if (method == "pois" ) mod$inits$theta =   log(ifelse(data$Y == 0, sum(data$Y) / sum(data$n), data$Y / data$n))
    initmiss = c(initmiss, "theta")
  }
  if (is.null(mod$inits$beta)) {
    if (method == "binom") mod$inits$beta = logit(apply(data$Y, 1:2, sum) / apply(data$n, 1:2, sum))
    if (method == "pois" ) mod$inits$beta =   log(apply(data$Y, 1:2, sum) / apply(data$n, 1:2, sum))
    initmiss = c(initmiss, "beta")
  }
  if (is.null(mod$inits$tau2)) {
    mod$inits$tau2 = matrix(0.01, nrow = Ng)
    initmiss = c(initmiss, "tau2")
  }
  if (is.null(mod$inits$G)) {
    mod$inits$G = diag(Ng) / 7
    initmiss = c(initmiss, "G")
  }
  if (is.null(mod$inits$Gt)) {
    mod$inits$Gt = array(diag(1, nrow = Ng), dim = c(Ng, Ng, Nt))
    initmiss = c(initmiss, "Gt")
  }
  if (is.null(mod$inits$z)) {
    mod$inits$z = mod$inits$theta - array(mod$inits$beta, dim = mod$params$dNd)
    for (g in 1:Ng) {
      for (t in 1:Nt) {
        mod$inits$z[g, t, ] = mod$inits$z[g, t, ] - mean(mod$inits$z[g, t, ])
      }
    }
    initmiss = c(initmiss, "z")
  }
  if (is.null(mod$inits$rho)) {
    mod$inits$rho = matrix(0.95, nrow = Ng)
    if (rho_up) initmiss = c(initmiss, "rho")
  }
  if (!is.null(initmiss)) cat("The following objects were created using defaults in 'inits':", paste(initmiss, collapse = " "), "\n")
  inicheck(mod, mod$inits)
  plot(1:Nt, (data$Y / data$n)[1, , which.max(data$n[1, 1, ])] * 1e5, xlab = "Time", ylab = "Rate", main = "Change in Rates Across Time")
  cat("Model set up!\n")
  mod
}
