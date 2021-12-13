#' Runs a multivariate point null test. This function returns
#' an approximate p-value for the specified test statistic.
#'
#' @param obs_data The observed data to be used for finding the optimal
#' norm (training), and finding the test statistic (testing).  Similar to
#' above, each row is an observation and each column corresponds to either
#' the outcome (first column) or a covariate.
#' @param param_est Function used to estimate the parameter and corresponding
#' influence curve.
#' @param control List used to define controls for test.
#' @return Either a p-value, or a p-value and other information, depending on
#' the value of control\$more_info.
#'
#' @export

mv_pn_test <- function(obs_data, param_est = NULL,
                    control = test.control()) {
  init_est <- calc_gam_star(
    obs_data = obs_data, param_est = param_est,
    control = control, lm_dst = NULL,
    return_lmd = TRUE
  )
  gam_star_n <- init_est$gam_star_n
  ic_est <- init_est$ic_est
  param_ses <- init_est$param_ses
  ts_ld_bs_samp <- control$ts_ld_bs_samp
  oth_ic_info <- init_est$oth_ic
  if (control$ld_est_meth == "par_boot") {
    sim_ts_mat <- matrix(stats::rnorm(ts_ld_bs_samp * nrow(ic_est)),
                         nrow = ts_ld_bs_samp)
    f_e_lm_dstr <- gen_boot_sample(sim_ts_mat, ic_est,
                                   center = TRUE, rate = "rootn")
    ts_lim_dist <- rep(NA, ts_ld_bs_samp)
    for (bs_idx in 1:ts_ld_bs_samp) {
      sub_data <- f_e_lm_dstr[bs_idx, , drop = FALSE]
      par_boot_est <- calc_gam_star(
        obs_data = sub_data,
        param_est = function(x, ...)list("est" = apply(x, 2, mean)),
        control = control, lm_dst = init_est$lm_dst,
        return_lmd = FALSE)
      ts_lim_dist[bs_idx] <- par_boot_est$gam_star_n
    }
  }else if (control$ld_est_meth == "perm") {
    ts_lim_dist <- rep(NA, ts_ld_bs_samp)
    num_obs <- nrow(obs_data)
    for (perm_idx in seq(ts_ld_bs_samp)) {
      y_idx <- sample(seq(num_obs), replace = FALSE)
      perm_data <- obs_data
      perm_data[, 1] <- perm_data[y_idx, 1]
      perm_est <- calc_gam_star(obs_data = perm_data, param_est = param_est,
                            control = control, lm_dst = NULL,
                            return_lmd = FALSE)
      ts_lim_dist[perm_idx] <- perm_est$gam_star_n
    }
  }
  p_val <- mean(as.integer(gam_star_n > ts_lim_dist))
  if (control$show_hist) {
    browser()
    graphics::hist(ts_lim_dist)
    graphics::abline(v = gam_star_n, col = "red")
  }
  if (!is.null(control$more_info)) {
    chsn_tbl <- vapply(control$pos_lp_norms,
                       function(x) mean(x == init_est$chsn_norms),
                       FUN.VALUE = -99)
    var_mat <- t(ic_est) %*% ic_est / nrow(ic_est)
    return(list("pvalue" = p_val,
                "test_stat" = gam_star_n,
                "test_st_eld" = ts_lim_dist,
                "chosen_norm" = chsn_tbl,
                "param_ests" = as.vector(init_est$param_ests),
                "param_ses" = param_ses,
                "var_mat" = var_mat,
                "oth_ic_inf" = oth_ic_info
                ))
  }else{
    return(c("p-value" = p_val))
  }
}
