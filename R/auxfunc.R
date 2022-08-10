getislands = function(mod) {
  check = 1:mod$params$dNd[3]
  network_true = list()
  p = 0
  while(length(check)) {
    network = NULL
    net     = check[1]
    network = c(net, mod$nb$neigh[[net]])
    if (!any(mod$nb$neigh[[net]] == 0)) {
      i = 2
      while(i <= length(network)) {
        network = append(network, mod$nb$neigh[[network[i]]][!(mod$nb$neigh[[network[i]]] %in% network)])
        net = append(net, network[i])
        i = i + 1
      }
    }
    p = p + 1
    network_true[[p]] = sort(net)
    check = check[-which(check %in% net)]
  }
  return(network_true)
}

nbcheck = function(mod) {
  # Check that all neighbor list objects are present
  chk  = c("num", "neigh")
  miss = sapply(1:length(chk), \(x) !any(names(mod$nb) == chk[x]))
  if (sum(miss)) stop("One or more objects missing from list 'nb': ", paste(chk[miss], collapse = ", "))

  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'nb'
  chk = which(!(names(mod$nb) %in% c("num", "neigh")))
  if (length(chk)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'nb':", paste(names(mod$nb)[chk], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) warning(paste(warnct, "warning(s) found in list 'nb':\n", paste(warnout, collapse = "\n ")))

  # Check for errors
  errout = NULL
  errct  = 0

  # Counties have no neighbors
  no_nb = which(sapply(seq_along(mod$nb$neigh), \(x) any(mod$nb$neigh[[x]] == 0)))
  if (length(no_nb)) {
    errct  = errct + 1
    errtxt = paste(errct, ": the following counties have no neighbors:", paste0(paste(no_nb, collapse = " "), "."), "Refer to documentation for fixes\n")
    errout = c(errout, errtxt)
  }
  # Neighbor data doesn't match shape data
  if (!is.null(mod$shp)) {
    shp_length = ifelse(is.null(nrow(mod$shp)), length(mod$shp), nrow(mod$shp))
    if (shp_length != length(mod$nb$neigh) | shp_length != length(mod$nb$num)) {
      errct  = errct + 1
      errtxt = paste(errct, ": Length of neighbor data does not match length of shape data. Ensure that the length/number of rows of input shapefile match the length of neighbor data\n")
      errout = c(errout, errtxt)
    }
  } else {
    if (dim(mod$data$Y)[3] != length(mod$nb$neigh) | dim(mod$data$Y)[3] != length(mod$nb$num)) {
      errct  = errct + 1
      errtxt = paste(errct, ": Length of neighbor data does not match length of data. Ensure that the length of event/population data match the length of neighbor data\n")
      errout = c(errout, errtxt)
    }
  }
  if (errct) stop(paste(errct, "error(s) found in list 'nb':\n", paste(errout, collapse = "\n ")))
}

datcheck = function(mod) {
  # Check that all data is present
  chk = c("Y", "n")
  miss = sapply(1:length(chk), \(x) !any(names(mod$data) == chk[x]))
  if (sum(miss)) stop("One or more objects missing from list 'data': ", paste(chk[miss], collapse = ", "))

  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'data'
  chk = which(!(names(mod$data) %in% c("Y", "n", "shp")))
  if (length(chk)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'data':", paste(names(mod$data)[chk], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) warning(paste(warnct, "warning(s) found in list 'data':\n", paste(warnout, collapse = "\n ")))

  # Check for errors
  errout = NULL
  errct  = 0
  # Dimensions of Y and n are not the same
  if (any(dim(mod$data$Y) != dim(mod$data$n))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Data not same dimensions. Ensure dim(Y) == dim(n)")
    errout = c(errout, errtxt)
  }
  if (errct) stop(paste(errct, "error(s) found in list 'data':\n", paste(errout, collapse = "\n ")))
}

pricheck = function(mod) {
  Ng = mod$params$dNd[1]; Nt = mod$params$dNd[2]
  # Check that all priors are present
  chk = c("eta", "Sig_b_i", "G0_i", "at", "bt", "nu", "nu_0", "gamma")
  if (mod$params$rho_up) chk = c(chk, "a_rho", "b_rho", "delta")
  miss = sapply(1:length(chk), \(x) !any(names(mod$priors) == chk[x]))
  if (sum(miss)) stop("One or more objects missing from list 'priors': ", paste(chk[miss], collapse = ", "))

  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'priors'
  chk = which(!(names(mod$priors) %in% c("eta", "Sig_b_i", "G0_i", "at", "bt", "nu", "nu_0", "gamma")))
  if (mod$params$rho_up) chk = which(!(names(mod$priors) %in% c("eta", "Sig_b_i", "G0_i", "at", "bt", "nu", "nu_0", "gamma",  "a_rho", "b_rho", "delta")))
  if (length(chk)) {
    warnct = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'priors':", paste(names(mod$priors)[chk], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) warning(paste(warnct, "warning(s) found in list 'priors':\n", paste(warnout, collapse = "\n ")))

  # Check for errors
  errout = NULL
  errct  = 0
  # Length of eta is not Ng * Nt
  if (length(mod$priors$eta) != Ng * Nt) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'eta' is invalid length. Ensure length(eta) == Ng * Nt or use default value")
    errout = c(errout, errtxt)
  }
  # Sig_b_i is not Ng * Nt x Ng * Nt
  if (any(dim(mod$priors$Sig_b_i) != rep(Ng * Nt, 2))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'Sig_b_i' is not a p x p matrix. Ensure dim(Sig_b_i) == Ng * Nt x Ng * Nt or use default value")
    errout = c(errout, errtxt)
    # Sig_b_i is Ng * Nt x Ng * Nt but not positive definite
  } else if (!isSymmetric(mod$priors$Sig_b_i)) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'Sig_b_i' is not symmetric. Ensure Sig_b_i is symmetric or use default value")
    errout = c(errout, errtxt)
  } else if (any(round(eigen(mod$priors$Sig_b_i)$values, 10) <= 0)) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'Sig_b_i' is not positive definite. Ensure Sig_b_i is positive definite or use default value")
    errout = c(errout, errtxt)
  }
  # nu is less than (Ng - 1)
  if (mod$priors$nu < mod$params$dNd[1] - 1) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'nu' is too small. Ensure nu >= Ng - 1 or use default value")
    errout = c(errout, errtxt)
  }
  # nu_0 is less than (Ng - 1)
  if (mod$priors$nu_0 < mod$params$dNd[1] - 1) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'nu_0' is too small. Ensure nu_0 >= Ng - 1 or use default value")
    errout = c(errout, errtxt)
  }
  # at is not positive
  if (mod$priors$at <= 0) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'at' is not positive. Ensure at > 0 or use default value")
    errout = c(errout, errtxt)
  }
  # bt is not positive
  if (mod$priors$bt <= 0) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'bt' is not positive. Ensure bt > 0 or use default value")
    errout = c(errout, errtxt)
  }
  # gamma is not correct dim
  if (any(dim(mod$priors$gamma) != dim(mod$data$Y))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'gamma' is invalid dimensions Ensure dim(gamma) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # any element of gamma is not positive
  if (any(mod$priors$gamma <= 0)) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'gamma' has non-positive values. Ensure all values of gamma are positive or use default value")
    errout = c(errout, errtxt)
  }
  if (mod$params$rho_up) {
    # any element of delta is not positive
    if (any(mod$priors$delta <= 0)) {
      errct  = errct + 1
      errtxt = paste(errct, ": Param 'delta' has non-positive values. Ensure all values of delta are positive or use default value")
      errout = c(errout, errtxt)
    }
    # a_rho is not positive
    if (mod$priors$a_rho <= 0) {
      errct  = errct + 1
      errtxt = paste(errct, ": Param 'a_rho' is not positive. Ensure a_rho > 0 or use default value")
      errout = c(errout, errtxt)
    }
    # b_rho is not positive
    if (mod$priors$b_rho <= 0) {
      errct  = errct + 1
      errtxt = paste(errct, ": Param 'b_rho' is not positive. Ensure b_rho > 0 or use default value")
      errout = c(errout, errtxt)
    }
  }
  if (errct) stop(paste(errct, "error(s) found in list 'priors':\n", paste(errout, collapse = "\n ")))
}

inicheck = function(mod, inits) {
  Ng = mod$params$dNd[1]; Nt = mod$params$dNd[2]
  # Check that all inits are present
  chk  = c("theta", "beta", "tau2", "G", "Gt", "z", "rho")
  miss = sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))

  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'inits'
  chk = which(!(names(inits) %in% c("theta", "beta", "tau2", "G", "Gt", "z", "rho")))
  if (length(chk)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'inits':", paste(names(inits)[chk], collapse = ", "))
    warnout = c(warnout, warntxt)
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # beta is not Ng x Nt
  if (any(dim(beta) != c(Ng, Nt))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'beta' is invalid dimensions. Ensure dim(beta) == Ng x Nt or use default value")
    errout = c(errout, errtxt)
  }
  # theta is not same dimensions as Y
  if (any(dim(inits$theta) != dim(mod$data$Y))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'theta' is invalid dimensions. Ensure dim(theta) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # z is not same dimensions as Y
  if (any(dim(inits$z) != dim(mod$data$Y))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Param 'z' is invalid dimensions. Ensure dim(z) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  if (errct) stop(paste(errct, "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
}

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

init = function(
	data,
	nb = NULL,
	shp = NULL,
	priors = NULL,
	inits = NULL,
	name = "output",
	dir = getwd(),
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
  mod$params$its  = 0
  mod$params$dNd  = dim(mod$data$Y)
  mod$params$I    = length(getislands(mod))
  mod$params$name = name
  mod$params$dir  = dir
  mod$nb$neigh    = sapply(mod$nb$neigh, \(x) x - 1) # translate neighbor information to zero-index
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
    if (method == "binom") mod$inits$theta = .logit(ifelse(data$Y == 0 | data$n == 0 | data$Y >= data$n, sum(data$Y) / sum(data$n), data$Y / data$n))
    if (method == "pois" ) mod$inits$theta =    log(ifelse(data$Y == 0, sum(data$Y) / sum(data$n), data$Y / data$n))
    initmiss = c(initmiss, "theta")
  }
  if (is.null(mod$inits$beta)) {
    if (method == "binom") mod$inits$beta = .logit(apply(data$Y, 1:2, sum) / apply(data$n, 1:2, sum))
    if (method == "pois" ) mod$inits$beta =    log(apply(data$Y, 1:2, sum) / apply(data$n, 1:2, sum))
    initmiss = c(initmiss, "beta")
  }
  if (is.null(mod$inits$tau2)) {
    mod$inits$tau2 = matrix(0.01, nrow = Ng)
    initmiss   = c(initmiss, "tau2")
  }
  if (is.null(mod$inits$G)) {
    mod$inits$G  = diag(Ng) / 7
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
  #plot(1:Nt, (data$Y / data$n)[1, , which.max(data$n[1, 1, ])] * 1e5, xlab = "Time", ylab = "Rate", main = "Change in Rates Across Time")
  cat("Model set up!\n")
  return(mod)
}

.logit = function(x) return(log(x / (1 - x)))
.expit = function(x) return(1 / (1 + exp(-x)))

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

getsuff = function(mod, param, burn) {
  dirname = paste0(getwd(), "/", mod$params$name)
  expr = paste0("^", mod$params$name, "_", param, "_[1-9]+[0-9]*\\.txt$")
  file_suff = list.files(paste0(dirname, "/", param))
  file_suff = sub(paste0("^", mod$params$name, "_.+_"), "", file_suff)
  file_suff = sub(".txt$", "", file_suff)
  file_suff = unique(sort(as.numeric(file_suff)))
  file_suff = file_suff[file_suff > burn]
  return(file_suff)
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
load_model = function(name, dir = getwd()) return(readRDS(paste0(dir, "/", name, "/mod_", name, ".Rds")))

#' Load Samples From Storage
#'
#' @param mod a model object created with the \code{init()} function.
#' @param params a string vector of parameters to load samples for. \code{all} specifies all parameters.
#' @param thin amount of thinning to be done on samples. Loads in every \code{thin} iteration of the samples.
#' Must multiply evenly into 100.
#' @param burn amount of burn-in to be done on samples. Discards first \code{burn} iterations of the samples.
#' @return A list containing samples for each parameter of interest.
#' @examples
#' \dontrun{
#' load_samples(mod_nc)
#' load_samples(mod_nc, params = "all")
#' load_samples(mod = mod_nc, params = "theta", thin = 10, burn = 2e3)
#' load_samples(mod_nc, c("theta", "tau2", "Gt"), burn = 2e3)
#' }
load_samples = function(mod, params = c("all", names(mod$inits)), thin = 1, burn = 0) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all" %in% params)  params = names(mod$inits)
  if (!mod$params$rho_up) params = params[which(params != "rho")]
  if (burn >= mod$params$its) stop("'burn' too big. Choose a 'burn' that is smaller than the total number of samples")
  if (100  %% thin != 0) stop("'thin' must multiply evenly into 100. Please choose another value for 'thin'")
  if (burn %% 100  != 0) stop("'burn' must be a multiple of 100. Please choose another value for 'burn'")

  output = .load_samples(mod, params, burn, thin, getsuff(mod, params[1], burn))
  output = lapply(output, simplify2array)

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
  return(output)
}

# add documentation for acceptance ratio stuff
acceptance_ratio = function(mod, params = c("all", "theta", "rho"), burn = 0) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all" %in% params) params = c("theta", "rho")
  if (!mod$params$rho_up) params = params[which(params != "rho")]
  accept = NULL
  if ("theta" %in% params) accept$theta = acceptance_ratio_cube(mod, getsuff(mod, "theta", burn), burn)
  if ("rho"   %in% params) accept$rho   = acceptance_ratio_mat (mod, getsuff(mod, "rho"  , burn), burn)
  return(accept)
}


# Maybe consoldiate this into recover_samples function???

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
get_medians = function(rec_samples, params = c("all", names(rec_samples)), ci = 0.95) {
  params = match.arg(params, several.ok = TRUE)
  params = unique(params)
  if ("all" %in% params) params = names(rec_samples)
  est = NULL
  for (nam in params) {
    d = 1:(length(dim(rec_samples[[nam]])) - 1)
    est[[nam]] = apply(rec_samples[[nam]], d, median)
  }
  return(est)
}
