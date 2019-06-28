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

mag_for_pow <- function(boot_data, dir, lp = 2, power = 0.8, nf_quant){
  n_covs <- ncol(boot_data)
  n_obs  <- nrow(boot_data)
  n_dir <- dir/(sqrt(sum(dir ** 2)))
  roots <- rep(NA, n_obs)
  if(lp == "max"){
    prd_vals <- sweep(boot_data, 2, n_dir, FUN = "*")
    maxes <- apply(prd_vals, 1, abs_max)
    smallish <- quantile(maxes, 1 - power)
    return(nf_quant/smallish)
  }else{
    lp <- as.numeric(lp)
  }
  poly <- rep(NA, lp + 1)
  dir_mat <- matrix(NA, nrow = lp + 1, ncol = n_covs)
  dir_mat[1, ] <- rep(1, n_covs)
  for(i in 2:(lp + 1)){
    dir_mat[i, ] <- dir_mat[i - 1, ] * n_dir
  }
  for(root_idx in 1:n_obs){
    obs <- boot_data[root_idx, ]
    for(k_idx in 0:lp){
      poly[k_idx + 1] <- choose(lp, k_idx) * sum(dir_mat[k_idx + 1, ] * (obs ** (lp - k_idx)))
    }
    poly[1] <- poly[1] - nf_quant ** lp
     all_rts <- polyroot(poly)
    # print(min(abs(Im(all_rts))))
     pos_root <- which(Re(all_rts) > 0)
     apr <- all_rts[pos_root]
     real_root <- which(abs(Im(apr)) < 0.00001)
     if(length(real_root) != 1){
       roots[root_idx] <- -1
     }else{
       roots[root_idx] <- Re(apr[real_root])
     }
  }
  return(quantile(roots, power))
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









