#' A helper function that calculates the estimated p-value for a given
#' observed alternative and a given norm.
#' @param mc_limit_dstr the simulated data draw from
#' the limiting distribution under the null
#' @param dir the observed estimate of the parameter
#' @param norms_idx the index for the norm used
#' @param norm_type the type of norm used
#' @param ... additional arguments that may be passed to
#' pval_for_mag, but which will be ignored.
#'
#' @return The p-value for a test with estimate given by
#' \code{dir} and limiting distribution given by \code{mc_limit_dstr}.
#'
#' @export

pval_for_mag <- function(mc_limit_dstr, dir, norms_idx = 2,
                         norm_type = "lp", ...) {
  num_norms <- length(norms_idx)
  est_pvals <- rep(NA, num_norms)
  distr <- matrix(NA, nrow = nrow(mc_limit_dstr), ncol = num_norms)
  if (norm_type == "ssq") {
    num_obs <- nrow(mc_limit_dstr)
    trans_dir <- cumsum(dir ** 2)
    for (obs_idx in 1:num_obs) {
      distr[obs_idx, ] <- cumsum(sort(mc_limit_dstr[obs_idx, ] ** 2,
                                      decreasing = TRUE))
    }
    for (lp_idx in 1:num_norms) {
      trans_est <- trans_dir[norms_idx[lp_idx]]
      est_pvals[lp_idx] <- mean(as.numeric(
        distr[, norms_idx[lp_idx]] >= trans_est
        ))
    }
    return(est_pvals)
  }
  if (norm_type == "lp") {
    nrmd_distr <- apply(mc_limit_dstr, 1, l_p_norm, p = norms_idx[lp_idx], type = "lp")
    for (lp_idx in 1:num_norms) {
      trans_dir <- l_p_norm(dir, p = norms_idx[lp_idx], type = "lp")
      est_pvals[lp_idx] <- mean(as.numeric(nrmd_distr >= trans_dir))
    }
    names(est_pvals) <- norms_idx
    return(est_pvals)
  }
}
