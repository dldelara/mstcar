#' @export
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
#' @export
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
#' @export
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
#' @export
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
    warnct  = warnct + 1
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
#' @export
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

#' @export
getsuff = function(mod, param, burn) {
  dirname   = paste0(getwd(), "/", mod$params$name)
  expr      = paste0("^", mod$params$name, "_", param, "_[1-9]+[0-9]*\\.txt$")
  file_suff = list.files(paste0(dirname, "/", param))
  file_suff = sub(paste0("^", mod$params$name, "_.+_"), "", file_suff)
  file_suff = sub(".txt$", "", file_suff)
  file_suff = unique(sort(as.numeric(file_suff)))
  file_suff = file_suff[file_suff > burn]
  return(file_suff)
}
