#' Estimate the local acceptance rate
#'
#' This a helper function used to estimate the acceptance rate for a
#' simple norm based test under a sequence of local alternatives in
#' a given direction given the estimated limiting distribution.
#'
#' @param mc_limit_dstr MC draws from an estimate of the
#' (centered) limiting distribution. Data where columns correspond
#' to different covariates, and rows are independent observations.
#' @param norms_idx the index of the norms to be used (ideally integers).
#' @param dir a vector in the direction for which we wish to estimate power.
#' @param null_quants the cutoff values for the distribution under the null
#' for each lp norm.
#' @param norm_type string indicating the class of norms to select over
#' (sum of squares ("ssq") or lp norms ("lp") are currently supported).
#'
#' @return The estimated acceptance rate for the specified norm-based test
#' for a given local alternative.

accept_rate <- function(mc_limit_dstr, dir,
                        null_quants, norms_idx = 2,
                        norm_type = "lp") {
  num_norms <- length(norms_idx)
  est_acc <- rep(NA, num_norms)
  shift_distr <- sweep(mc_limit_dstr, 2, dir, "+")
  if (norm_type == "ssq") {
    num_obs <- nrow(mc_limit_dstr)
    for (obs_idx in 1:num_obs) {
      shift_distr[obs_idx, ] <- cumsum(sort(shift_distr[obs_idx, ] ** 2,
                                             decreasing = TRUE))
    }
    for (lp_idx in 1:num_norms) {
      est_acc[lp_idx] <- mean(as.numeric(shift_distr[, norms_idx[lp_idx]] <
                                            null_quants[lp_idx]))
    }
  }
  if (norm_type == "lp") {
    for (lp_idx in 1:num_norms) {
      norm_shift_distr <- apply(shift_distr, 1, l_p_norm,
                                p = norms_idx[lp_idx], type = "lp")
      est_acc[lp_idx] <- mean(as.numeric(norm_shift_distr <
                                           null_quants[lp_idx]))
    }
  }
  names(est_acc) <- norms_idx
  return(est_acc)
}
