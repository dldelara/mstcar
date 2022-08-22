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
