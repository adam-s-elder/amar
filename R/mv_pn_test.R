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
#' @return The test will always return the following output:
#'
#' - \code{pvalue}: The approximate value of the test statistic
#'
#' - \code{test_stat}: The approximate value of the test statistic
#'
#' - \code{test_st_eld}: The approximate limiting distribution of the test statistic
#'  (with length equal to \code{ts_ld_bs_samp}).
#'
#'- \code{chosen_norm}: A vector indicating which norm was chosen by
#' the adaptive test
#'
#'- \code{param_ests}: The parameter estimate.
#'
#'- \code{param_ses}: An estimate of the standard error of
#'  each element of \code{param_ests}
#'
#'- \code{oth_ic_inf}: Any other information provided by the \code{param_est}
#'  function when calculating the IC and parameter estimates.
#' Additional information may be returned by specifying it in the
#' test.control function:
#'
#' - If \code{"var_est"} is contained in \code{other_output}, the test output
#' will contain will have \code{var_mat} returned which is the empirical
#' second moment of the IC (equal asymptotically to the variance estimator).
#'
#' - If \code{"obs_data"} is contained in the \code{other_output}, the test
#' output will return the data passed to the testing function.
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' ## NOTE: More monte-carlo samples should be taken are taken here.  This is
#' ## only done to lower computation time.
#' test <- mv_pn_test(data.frame(y = rnorm(100), x = rnorm(100)),
#'                    ic.pearson, test.control(n_peld_mc_samples = 20,
#'                                             ts_ld_bs_samp = 20))
#'

mv_pn_test <- function(obs_data, param_est = NULL,
                       control = test.control()) {
  .checkargs(obs_data, param_est, control)
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
  chsn_tbl <- vapply(control$pos_lp_norms,
                     function(x) mean(x == init_est$chsn_norms),
                     FUN.VALUE = -99)
  out <- list("pvalue" = p_val,
              "test_stat" = gam_star_n,
              "test_st_eld" = ts_lim_dist,
              "chosen_norm" = chsn_tbl,
              "param_ests" = as.vector(init_est$param_ests),
              "param_ses" = param_ses,
              "oth_ic_inf" = oth_ic_info
  )
  if ("var_est" %in% control$other_output) {
    out$var_mat <- t(ic_est) %*% ic_est / nrow(ic_est)
  }
  if ("data" %in% control$other_output) {
    out$obs_data <- obs_data
  }
  return(out)
}

.checkargs <- function(obs_data, param_est, control) {
  default_control_args <- names(amp::test.control())
  passed_controls <- names(control)
  mis_args <- setdiff(default_control_args, passed_controls)
  if (length(mis_args) > 0) stop(
    paste0("Some control arguments have not been provided: ",
           paste0("\n ", mis_args, " ", collapse = ""),
           "\n To avoid this error, consider using amp::test.control()")
  )
  non_stand_args <- setdiff(passed_controls, default_control_args)
  if (length(non_stand_args) > 0) {
    fun_def <- deparse(param_est)
    unused_args <- rep(TRUE, length(non_stand_args))
    for (i in seq_len(length(unused_args))) {
      unused_args[i] <- !any(grepl(non_stand_args[i], fun_def))
    }
    if (any(unused_args)) {
      ua_names <- non_stand_args[unused_args]
      warning(
        paste0("Some arguments appear not to be used by the test.  These are: ",
               paste0("\n ", ua_names, " ", collapse = ""),
               "\n To avoid this warning, explicitly call these arguments ",
               "in your param_est function (rather than just ",
               "passing all control arguments to another function.")
      )
    }
  }
  if (!(control$ld_est_meth %in% c("par_boot", "perm"))) {
    stop("The control argument ld_est_meth must be either 'par_boot' or 'perm'")
    }
}
