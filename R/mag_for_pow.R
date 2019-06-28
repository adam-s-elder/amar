#' Helper function for est_local_pow.  This function
#' is used to estimate the magnitude needed for a
#' certain direction to achive a given power.
#'
#' @param boot_data bootstrap estimates of the (centered) limiting distribution.
#' Data where columns correspond to different covariates, and rows are independent observations.
#' @param lp the lp norm to be used (ideally it is an integer)
#' @param dir the direction for which we wish to find the magnitude needed to achieve
#' the speceified power.
#' @param power power which the local alternative is to achieve
#' for a magnitude that will be solved.
#'
#' @return Magnitude for the specified lp norm for a given local alternative.
#'
#' @export

abs_max <- function(x) max(abs(x))

mag_for_pow <- function(boot_data, dir, lp_nrms = 2, power = 0.8, nf_quants, nrm_type = "lp"){
  n_covs <- ncol(boot_data)
  n_obs  <- nrow(boot_data)
  roots <- rep(NA, n_obs)
  n_norms <- length(lp_nrms)
  norm_res <- rep(NA, n_norms)
  for(norm_idx in 1:n_norms){
    lp <- lp_nrms[norm_idx]
    nf_quant <- nf_quants[norm_idx]
    if(lp == "max"){
      many_mags <- rep(NA, n_obs)
      for(obs_idx in 1:n_obs){
        xs <- boot_data[obs_idx]
        all_opts <- pmax((nf_quant - xs) / dir, (-nf_quant - xs )/ dir)
        many_mags[obs_idx] <- min(pmax(0, all_opts))
      }
      mfp <- quantile(many_mags, power)
      norm_res[norm_idx] <- mfp
    }else{
      lp <- as.numeric(lp)
      for(root_idx in 1:n_obs){
        sing_obs <- boot_data[root_idx, ]
        roots[root_idx] <- find_mag(one_obs = sing_obs, dir = dir, cutoff = nf_quant,
                                    nrm_idx = lp, nrm_type = nrm_type)
      }
      norm_res[norm_idx] <- quantile(roots, power)
    }
  }
  names(norm_res) <- lp_nrms
  return(norm_res)
}


# xx <- seq(from = -3, to = 3, length.out = 1000)
# yy <- rep(NA, 1000)
# zz <- rep(NA, 1000)
# for(bbs in 1:1000) {
#   yy[bbs] <- (sum(poly * (xx[bbs] ** c(0:lp)))) ** (1/lp)
#   zz[bbs] <- l_p_norm(obs + xx[bbs] * n_dir, p = lp)
# }
# plot(xx, yy, type = "l", ylim = c(0, max(yy)))
# abline(h = nf_quant, col = "red")
# plot(xx, zz, type = "l", ylim = c(0, max(zz)))
# abline(h = nf_quant, col = "red")
# browser()









