#' Generate an estimate of the limiting distribution of
#' the vector covariate specific parameter estimates
#'
#' @param data observed data where columns correspond to different covariates, and rows are independent observations.
#' @param param The parameter of interest for example spearman correlation.
#' @param n the number of monte carlo draws. Default is 1000.
#' @return An n by p matrix in which each column corresponds to the
#' bootstrap parameter estimate for the corresponding covariate under
#' the hypothesis that all parameter estimates are zero.
#' @export
#'

est_limit_distr <- function(data, param, n = 1000){
  num_obs <- nrow(data)
  ic_ests <- est_influence_pearson(data)
  eps_mat <- matrix(rnorm(n * num_obs), ncol = num_obs, nrow = n)
  cent_boot <- eps_mat %*% ic_ests / sqrt(num_obs)
  return(cent_boot)
}

