#' A function used to select which l_p norm to use for
#' a given test.
#'
#' @param e_lm_dstr is the estimated limiting distribution
#' @param pos_lp is a vector containing the possible
#' values of the lp norm to be considered
#' @param search_size how many random draws will be taken to
#' estimate the average power.
#' @param fxd_val which parameter (power or magnitude) you wish
#' to keep fixed when solvinf for the best lp norm.  If fxd_val
#' is "power", the criteria used to choose the lp norm will be
#' based on the magnitude of local alternative required to reach
#' a specific power.
#' @param compare_measure whether you wish to compare the average scinario
#' verus the worst case scenario.
#'
#' @export

select_lp <- function(e_lm_dstr, pos_lp, search_size = 1000, fxd_val, compare_measure){
  num_norms <- length(pos_lp)
  num_covs  <- ncol(e_lm_dstr)
  cutoffs   <- rep(NA, num_norms)
  if(fxd_val == "power"){
    mag_meas <- mag_for_pow ; worst <- max ; best <- min ; cns <- 0.8
  }else if (fxd_val == "mag"){
    mag_meas <- pow_for_mag ; worst <- min ; best <- max ; cns <- 2
  }
  for(bb in 1:num_norms){
    lp_norm_distr <- apply(e_lm_dstr, 1, l_p_norm, p = pos_lp[bb])
    cutoffs[bb]   <- stats::quantile(lp_norm_distr, 0.95)
      }
  if (compare_measure == "mean"){
    search_dirs <- r_unif_sphere(num_covs, search_size)
    res_mat <- matrix(NA, nrow = length(pos_lp), ncol = search_size)
    for(jj in 1:search_size){
      for(p_idx in 1:num_norms){
        res_mat[p_idx, jj] <- mag_meas(e_lm_dstr, dir = search_dirs[jj, ],
                                       lp = pos_lp[p_idx], cns, nf_quant = cutoffs[p_idx])
      }
    }
   perfs <- apply(res_mat, 1, mean)
   names(perfs) <- pos_lp
  }else if (compare_measure == "extr"){

  }else{
    warning("compare_measure must either be mean or extr")
    return(NULL)
  }
  best_norm_idx <- which(best(perfs) == perfs)
  return(pos_lp[best_norm_idx])
}


