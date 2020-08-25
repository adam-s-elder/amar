#' Helper function for est_local_pow.  This function
#' is used to estimate the power of a
#' certain direction given the magnitude of that direction.
#'
#' @param boot_data bootstrap estimates of the (centered) limiting distribution.
#' Data where columns correspond to different covariates, and rows are independent observations.
#' @param lp the lp norms to be used (ideally integers).
#' @param dir a vector in the direction for which we wish to estimate power.
#' @param nf_quants the cutoff values for the distribution under the null and each lp norm.
#' @param nrm_type the class of norms to select over.
#'
#' @return Magnitude for the specified lp norm for a given local alternative.
#'
#' @export

pow_for_mag <- function(boot_data, dir, lp = 2, nf_quants, nrm_type = "lp"){
  num_norms <- length(lp)
  est_pows <- rep(NA, num_norms)
  shift_distr <- sweep(boot_data, 2, dir, "+")
  if(nrm_type == "ordl2"){
    num_obs <- nrow(boot_data)
    for(obs_idx in 1:num_obs){
      shift_distr[obs_idx, ] <- cumsum(sort(shift_distr[obs_idx, ] ** 2, decreasing = TRUE))
    }
    for(lp_idx in 1:num_norms){
      est_pows[lp_idx] <- mean(as.numeric(shift_distr[, lp[lp_idx]] >= nf_quants[lp_idx]))
    }
    return(est_pows)
  }
  if(nrm_type == "lp"){
    for(lp_idx in 1:num_norms){
      norm_shift_distr <- apply(shift_distr, 1, l_p_norm, p = lp[lp_idx], type = "lp")
      est_pows[lp_idx] <- mean(as.numeric(norm_shift_distr >= nf_quants[lp_idx]))
    }
    return(est_pows)
  }
}
