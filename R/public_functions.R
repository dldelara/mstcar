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
#' Gibbs Sampler Updates
#'
#' @param mod a model object created with the \code{init()} function.
#' @param n_iter number of iterations to sample. Must be a multiple of 100 and be ≥ \code{r}.
#' @param r number of iterations to store in each batch. Must be a multiple of 100 and multiply evenly into \code{n_iter}.
#' @examples
#' \dontrun{
#' samples(mod_nc, 6000)
#' samples(mod = mod_nc, n_iter = 1e4, r = 500)
#' }
#' @export
samples = function(mod, n_iter, r = 100) {
  dirname = paste0(mod$params$dir, "/", mod$params$name)
  while ((r %% 100) != 0) {
    r = as.numeric(readline(paste0("'r' must be a multiple of 100. Your r is currently ", r, ". Please specify a new r: ")))
  }
  while (((n_iter %% r) != 0) | (n_iter < r)) {
    if (n_iter < r) {
      n_iter = as.numeric(readline(paste0("'n_iter' must be \u2265 r. Your n_iter is currently ", n_iter, " and r is currently ", r, ". Please specify a new n_iter: ")))
    }
    if ((n_iter %% 100) != 0) {
      n_iter = as.numeric(readline(paste0("'n_iter' must be a multiple of 100. Your n_iter is currently ", n_iter, ". Please specify a new n_iter: ")))
    }
  }
  while (any(grepl(paste0("^.*/", mod$params$name, "_.*_[0-9]*.txt$"), list.files(dirname, recursive = TRUE))) & (mod$params$its == 0)) {
    delete = menu(
      c("No, let me change the name", "Yes", "Nevermind, exit for now"),
      title = paste("There are already output files with the name", mod$params$name, "in this directory. Running this sampler will delete these files. Run anyway?")
    )
    if (delete == 1) mod$params$name = readline("What would you like your new name to be? ")
    if (delete == 2) {
      unlink(paste0(dirname, "/", list.files(dirname, recursive = TRUE)[grepl("*.txt$", list.files(dirname, recursive = TRUE))]))
    }
    if (delete == 3) return(invisible())
  }
  if (!file.exists(dirname)) dir.create(dirname)
  pars = names(mod$inits)
  if (!mod$params$rho_up) pars = pars[which(pars != "rho")]
  for (par in pars) {
    if (!file.exists(paste0(dirname, "/", par))) dir.create(paste0(dirname, "/", par))
  }
  if (mod$params$its == 0) saveRDS(mod$inits, paste0(dirname, "/", mod$params$name, "_0.Rds"), version = 2)
  l = n_iter / r
  for (n_loop in 1:l) {
    gibbs_sampler(mod, r, n_loop, l)
    saveRDS(mod, paste0(dirname, "/mod_", mod$params$name, ".Rds"), version = 2)
  }
}

#' Load Model From Directory
#'
#' @param name the name assigned to the model in the \code{init} function.
#' @param dir the directory where the model exists.
#' @return A list containing model data, neighbor data, hyperparameters,
#' and initial output
#' @examples
#' \dontrun{
#' load_model("heart_nc")
#' load_model(name = "heart_nc", dir = getwd())
#' }
#' @export
load_model = function(name, dir = getwd()) {
  mod = readRDS(paste0(dir, "/", name, "/mod_", name, ".Rds"))
  mod$params$dir = gsub('(/)\\1+', '\\1', dir)
  mod
}
#' Load Samples From Storage
#'
#' @param mod a model object created with the \code{init()} function.
#' @param burn amount of burn-in to be done on samples. Discards first \code{burn} iterations of the samples.
#' @param thin amount of thinning to be done on samples. Loads in every \code{thin} iteration of the samples.
#' Must multiply evenly into 100.
#' @param params a string vector of parameters to load samples for. \code{all} specifies all parameters.
#' @return A list containing samples for each parameter of interest.
#' @examples
#' \dontrun{
#' load_samples(mod_nc)
#' load_samples(mod_nc, params = "all")
#' load_samples(mod = mod_nc, params = "theta", thin = 10, burn = 2e3)
#' load_samples(mod_nc, c("theta", "tau2", "Gt"), burn = 2e3)
#' }
#' @export
load_samples = function(mod, burn = 0, thin = 1, params = c("all", names(mod$inits))) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all" %in% params)  params = names(mod$inits)
  if (!mod$params$rho_up) params = params[which(params != "rho")]
  if (burn >= mod$params$its) {
    stop("'burn' too big. Choose a 'burn' that is smaller than the total number of samples")
  }
  if (100  %% thin != 0) {
    stop("'thin' must multiply evenly into 100. Please choose another value for 'thin'")
  }
  if (burn %% 100  != 0) {
    stop("'burn' must be a multiple of 100. Please choose another value for 'burn'")
  }
  cat("About to run get_output() function...\n\n")
  output = get_output(mod, burn, thin, params, getsuff(mod, params[1], burn))
  cat("Finished running get_output()!\n\n")
  output = lapply(output, simplify2array)
  cat("Created output list with dimensions:\n")
  print(lapply(output, dim))

  if ("tau2" %in% params) output$tau2 = output$tau2[1, , ]
  if ("rho"  %in% params) output$rho  = output$rho [1, , ]
  if (!is.null(dimnames(mod$data$Y))) {
    dims = c(dimnames(mod$data$Y), list(seq(burn + thin, mod$params$its, by = thin)))
    if ("theta" %in% params) dimnames(output$theta) = dims
    if ("beta"  %in% params) dimnames(output$beta ) = dims[-3]
    if ("tau2"  %in% params) dimnames(output$tau2 ) = dims[c(1, 4)]
    if ("Gt"    %in% params) dimnames(output$Gt   ) = dims[c(1, 1, 2, 4)]
    if ("G"     %in% params) dimnames(output$G    ) = dims[c(1, 1, 4)]
    if ("z"     %in% params) dimnames(output$z    ) = dims
    if ("rho"   %in% params) dimnames(output$rho  ) = dims[c(1, 4)]
  }
  output
}
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
  if ("all" %in% params) params = names(rec_samples)
  est = NULL
  for (nam in params) {
    d = 1:(length(dim(rec_samples[[nam]])) - 1)
    est[[nam]] = apply(rec_samples[[nam]], d, quantile, c(0.5, (1 - ci) / 2, 1 - (1 - ci) / 2))
  }
  est
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
