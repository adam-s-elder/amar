#' This function is used to estimate the magnitude needed for a
#' certain direction to achieve 80\% power for a proposed alternative.
#'
#' @param mc_limit_dstr MC draws from the estimated (centered)
#' limiting distribution. Columns correspond to different covariates,
#' and rows are independent observations.
#' @param norms_idx the norms to be considered.
#' @param norm_type the class of norms to be used
#' @param dir the direction for which we wish to find the magnitude needed
#' to achieve 80\% power.
#' @param null_quants ninety five percent quantiles for each of the
#' different norms
#'
#' @return Magnitude for the specified lp norm for a given local alternative.
#'
#' @export

mag_for_pow <- function(mc_limit_dstr, dir, norms_idx = 2,
                        null_quants, norm_type = "lp") {
  n_obs  <- nrow(mc_limit_dstr)
  roots <- rep(NA, n_obs)
  n_norms <- length(norms_idx)
  norm_res <- rep(NA, n_norms)
  for (norm_idx in 1:n_norms) {
    lp <- norms_idx[norm_idx]
    nf_quant <- null_quants[norm_idx]
    if (lp == "max") {
      many_mags <- rep(NA, n_obs)
      for(obs_idx in 1:n_obs){
        xs <- mc_limit_dstr[obs_idx]
        all_opts <- pmax((nf_quant - xs) / dir, (-nf_quant - xs) / dir)
        many_mags[obs_idx] <- min(pmax(0, all_opts))
      }
      mfp <- stats::quantile(many_mags, 0.8)
      norm_res[norm_idx] <- mfp
    }else{
      lp <- as.numeric(lp)
      for (root_idx in 1:n_obs) {
        sing_obs <- mc_limit_dstr[root_idx, ]
        roots[root_idx] <- find_mag(one_obs = sing_obs, dir = dir,
                                    cutoff = nf_quant, nrm_idx = lp,
                                    nrm_type = norm_type)
      }
      norm_res[norm_idx] <- stats::quantile(roots, 0.8)
    }
  }
  names(norm_res) <- norms_idx
  return(norm_res)
}
