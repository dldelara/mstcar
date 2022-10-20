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
  mod$params$dir = path.expand(gsub('(/)\\1+', '\\1', dir))
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
  output = get_output(mod, burn, thin, params, getsuff(mod, params[1], burn))
  output = lapply(output, simplify2array)

  if ("tau2" %in% params) output$tau2 = t(output$tau2[1, , ])
  if ("rho"  %in% params) output$rho  = t(output$rho [1, , ])
  if (!is.null(dimnames(mod$data$Y))) {
    dims = c(dimnames(mod$data$Y), list(seq(burn + thin, mod$params$its, by = thin)))
    if ("theta" %in% params) dimnames(output$theta) = dims
    if ("beta"  %in% params) dimnames(output$beta ) = dims[-3]
    if ("tau2"  %in% params) dimnames(output$tau2 ) = dims[c(4, 1)]
    if ("Gt"    %in% params) dimnames(output$Gt   ) = dims[c(1, 1, 2, 4)]
    if ("G"     %in% params) dimnames(output$G    ) = dims[c(1, 1, 4)]
    if ("z"     %in% params) dimnames(output$z    ) = dims
    if ("rho"   %in% params) dimnames(output$rho  ) = dims[c(4, 1)]
  }
  output
}
getsuff = function(mod, param, burn) {
  dirname   = paste0(mod$params$dir, "/", mod$params$name)
  expr      = paste0("^", mod$params$name, "_", param, "_[1-9]+[0-9]*\\.txt$")
  file_suff = list.files(paste0(dirname, "/", param))
  file_suff = sub(paste0("^", mod$params$name, "_.+_"), "", file_suff)
  file_suff = sub(".txt$", "", file_suff)
  file_suff = unique(sort(as.numeric(file_suff)))
  file_suff = file_suff[file_suff > burn]
  return(file_suff)
}
