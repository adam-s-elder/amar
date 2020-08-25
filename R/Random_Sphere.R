#' This function is used to generate random observations from the
#' unit ball (mainly used to generate random directions).
#' @param n_covs Number of covariates
#' @param n_obs Number of observations
#'
#' @export


r_unif_sphere <- function(n_covs, n_obs){
  obs <- matrix(stats::rnorm(n_covs * n_obs), nrow = n_covs)
  return(t(apply(obs, 2, function(x) x/sqrt(sum(x ** 2)) )))
}
