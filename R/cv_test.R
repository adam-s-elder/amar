#' CV-test runs a cross-validation, parametric bootstrap test. This function returns
#' an approximate p-value for the specified test statistic.
#'
#' @param obs_data The observed data to be used for finding the optimal
#' norm (training), and finding the test statistic (testing).  Similar to
#' above, each row is an observation and each column corresponds to either
#' the outcome (first column) or a covariate.
#' @param param_est Function used to estimate the parameter and corresponding
#' influence curve.
#' @param control List used to define controls for test.
#' @return learned test statistic for a single fold of data
#'
#' @export

cv_test <- function(obs_data, param_est = NULL,
                    control = test.control()){
  test_type <- control$test_type
  if (control$nrm_type == "bonf") {test_type <- "bonf"}
  init_cv_est <- est_cv(obs_data = obs_data, param_est = param_est,
                        control = control, est_lm_distr = NULL,
                        return_lmd = TRUE, view_IC = control$show_hist)
  e_lm_dstr <- init_cv_est$est_lm_dstr
  t_s_f <- init_cv_est$t_s_f
  cv_est <- init_cv_est$cv_est
  ic_est <- init_cv_est$ic_est
  param_ests <- init_cv_est$param_ests
  param_sds <- init_cv_est$param_sds
  one_step_cor <- init_cv_est$correction
  num_folds <- control$num_folds
  ts_ld_bs_samp <- control$ts_ld_bs_samp
  if (test_type == "par_boot") {
    if (control$nrmlize) {
      f_e_lm_dstr <- matrix(
        stats::rnorm(num_folds * ts_ld_bs_samp * ncol(ic_est)),
        nrow = num_folds * ts_ld_bs_samp)
    }else{
      sim_ts_mat <- matrix(
        stats::rnorm(num_folds * ts_ld_bs_samp * nrow(ic_est)),
        nrow = num_folds * ts_ld_bs_samp)
      f_e_lm_dstr <- gen_boot_sample(sim_ts_mat, ic_est,
                                     center = TRUE, rate = "rootn")
    }
    ts_lim_dist <- rep(NA, ts_ld_bs_samp)
    for (bs_idx in 1:ts_ld_bs_samp) {
      sub_data <- f_e_lm_dstr[
        (num_folds * (bs_idx - 1) + 1):(num_folds * bs_idx), , drop = FALSE
        ]
      par_boot_cv_est <- est_cv(
        obs_data = sub_data,
        param_est = function(x)list("est" = apply(x, 2, mean)),
        control = control, est_lm_distr = e_lm_dstr,
        return_lmd = FALSE)
      ts_lim_dist[bs_idx] <- par_boot_cv_est$cv_est
    }
  }else if (test_type == "perm") {
    ts_lim_dist <- rep(NA, ts_ld_bs_samp)
    num_obs <- nrow(obs_data)
    for (perm_idx in seq(ts_ld_bs_samp)) {
      y_idx <- sample(seq(num_obs), replace = FALSE)
      perm_data <- obs_data
      perm_data[, 1] <- perm_data[y_idx, 1]
      perm_cv_est <- est_cv(obs_data = perm_data, param_est = param_est,
                            control = control, est_lm_distr = NULL,
                            return_lmd = FALSE)
      ts_lim_dist[perm_idx] <- perm_cv_est$cv_est
    }
  }else{
    ts_lim_dist <- NULL
  }

  if (test_type == "bonf") {
    zvals <- abs(param_ests) / param_sds
    p_val <- 2 * (1 - stats::pnorm(max(zvals))) * length(zvals)
  }else{
    init_p_val <- mean(as.integer(cv_est <= ts_lim_dist))
    if (control$test_stat_func %in% c("mag", "pval")) {
      p_val <- 1 - init_p_val
    }else{
      p_val <- init_p_val
    }
  }

  if (control$show_hist & p_val < 0.05) {
    browser()
    graphics::hist(ts_lim_dist)
    graphics::abline(v = cv_est, col = "red")
  }

  if (!is.null(control$more_info)) {
    chsn_tbl <- vapply(control$pos_lp_norms,
                       function(x) mean(x == init_cv_est$chsn_norms),
                       FUN.VALUE = -99)
    if (TRUE) { #more_info == "all"
      var_mat <- t(ic_est) %*% ic_est / nrow(ic_est)
    }else{
      ## Currently waiting to see if this adds too much
      ## to the list.  May remove
      var_mat <- NULL
    }

    return(list("pvalue" = p_val,
                "test_stat" = cv_est,
                "test_st_eld" = ts_lim_dist,
                "chosen_norm" = chsn_tbl,
                "param_ests" = param_ests,
                "onestep_cor" = one_step_cor,
                "param_sds" = param_sds,
                "var_mat" = var_mat
                ))
  }else{
    return(c("p-value" = p_val))
  }
}

