#' A helper function for \code{mv_pn_test}, calculating the test statistic
#' for both the vector of parameter estimates, and the draws from
#' the corresponding estimated limiting distribution.
#'
#' @param obs_data The observed data used to calculate the test statistic.
#' Each row is an observation and each column corresponds to either the
#' outcome (first column) or a covariate.
#' @param param_est A function used to estimate both they parameter of interest
#' and the IC of the corresponding estimator.
#' @param control A list providing control arguments for the function.
#' @param lm_dst A list containing an estimate of the limiting distribution and
#' corresponding norm specific test cutoffs if it is provided.
#' @param return_lmd Boolean for whether to return the estimated
#' limiting distribution.
#' @return Calculated test statistic for the given data.
#' @importFrom expm sqrtm
#'
#' @export

calc_gam_star <- function(obs_data, param_est, control,
                          lm_dst = NULL, return_lmd = FALSE) {
  pos_lp_norms <- control$pos_lp_norms
  if (!is.null(lm_dst)) {
    psi_est <- param_est(obs_data)$est
  }else{
    est_and_ic <- param_est(obs_data, what = "both", control = control)
    ic_est <- est_and_ic$ic
    # if (control$shrink_cov) { #Experimental; using a shrinkage estimator for IC:
    #   old_ic <- ic_est
    #   invisible(capture.output(
    #     ic_est <- nlshrinkadam::nlshrink_cov(ic_est, ret_ic = TRUE)
    #   ))
    #   ## The following line is added to make sure the new covariance estimate
    #   ## is scaled the same as the old IC.
    #   ic_est <- ic_est * sqrt(nrow(ic_est))
    # }
    if (control$standardize) {
        vcov_mat <- t(ic_est) %*% ic_est / nrow(obs_data)
        # print(vcov_mat)
        if (nrow(obs_data) < nrow(vcov_mat)) {
          stop("Support for standardization when p > n has not been added")
          # s_mat <- MASS::ginv(vcov_mat)
        }else{
          s_mat <- solve(vcov_mat)
        }
        # print(s_mat)
        sqrt_mat <- Re(expm::sqrtm(s_mat))
        psi_est <- sqrt_mat %*% as.vector(est_and_ic$est)
        ic_est <- diag(length(psi_est)) * sqrt(nrow(ic_est))
        zero_mat <- matrix(0, nrow = nrow(obs_data) - nrow(ic_est),
                           ncol = ncol(ic_est))
        ic_est <- rbind(ic_est, zero_mat)
      } else {
        psi_est <- est_and_ic$est
    }
    n_peld_mc_samples <- control$n_peld_mc_samples
    norm_mat  <- matrix(
      stats::rnorm(n_peld_mc_samples * nrow(ic_est)), nrow = n_peld_mc_samples
      )
    e_lm_dstr <- gen_boot_sample(
      norm_mat, ic_est, center = TRUE, rate = "rootn"
    )
    num_norms <- length(pos_lp_norms)
    cutoff_vals <- rep(NA, num_norms)
      for (nrm_idx in 1:num_norms) {
        normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                                p = pos_lp_norms[nrm_idx],
                                type = control$nrm_type)
        if (sum(is.na(normalized_obs)) != 0) browser()
        cutoff_vals[nrm_idx] <- stats::quantile(normalized_obs, 0.95)
      }
    lm_dst <- list(distr = e_lm_dstr, cutoffs = cutoff_vals)
  }
  par_est <- sqrt(nrow(obs_data)) * as.vector(psi_est)
  gammas <- control$perf_meas(
    mc_limit_dstr = lm_dst$distr, dir = par_est,
    null_quants = lm_dst$cutoffs,
    norms_idx = pos_lp_norms, norm_type = control$nrm_type
  )
  gam_star_n <- min(gammas)
  chsn_norm <- names(gammas)[which.min(gammas)]
  if (return_lmd) {
    param_ses <- apply(
      ic_est, 2, FUN = function(x) sqrt(mean((x - mean(x)) ** 2) / length(x))
      )
    # if (control$shrink_cov) {
    #   param_ses <-  apply(
    #     ic_est, 2, FUN = function(x) sqrt(mean((x) ** 2) / length(x))
    #   )
    # }
    if (exists("est_and_ic")) {
      oth_ic_inf <- est_and_ic[setdiff(names(est_and_ic), c("est", "ic"))]
      # oth_ic_inf$old_ic <- old_ic
    }else{
      oth_ic_inf <- NULL
    }
    return(list("gam_star_n" = gam_star_n, "chsn_norms" = chsn_norm,
                "lm_dst" = lm_dst, "t_s_f" = control$perf_meas,
                "ic_est" = ic_est, "param_ests" = psi_est,
                "param_ses" = param_ses, "oth_ic" = oth_ic_inf))
  }else{
    return(list("gam_star_n" = gam_star_n, "chsn_norms" = chsn_norm))
  }
}
