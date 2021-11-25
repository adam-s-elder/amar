#' This is an intermediate function inbetween the cv_test and one_fold_est functions which
#' allows for specification of the measure of performance.
#'
#' @param obs_data The observed data to be used for finding the optimal norm (training),
#' and finding the test statistic (testing).  Similar to above, each row is an observation and each
#' column corresponds to either the outcome (first column) or a covariate.
#' @param param_est A function used to estimate both they parameter of interest
#' and the IC of the corresponding estimator.
#' @param control A list providing control arguments for the function.
#' @param est_lm_distr An estimate of the limiting distribution
#' if it is provided.
#' @param return_lmd Boolean for whether to return the estimated
#' limiting distribution.
#' @param view_IC Boolean for whether to see summaries of the IC
#' @return learned test statistic for a single fold of data
#'
#' @export

est_cv <- function(obs_data, param_est, control, est_lm_distr = NULL,
                   return_lmd = FALSE, view_IC = FALSE) {
  n_bs_smp <- control$n_bs_smp
  perf_meas <- control$perf_meas
  if (!is.null(est_lm_distr)) {
    e_lm_dstr <- est_lm_distr
  }else{
    est_and_ic <- param_est(obs_data, what = "both", control = control)
    ic_est <- est_and_ic$ic
    psi_est <- est_and_ic$est
    if (control$nrmlize) {
      e_lm_dstr <- matrix(
        stats::rnorm(n_bs_smp * ncol(ic_est)), nrow = n_bs_smp
        )
    }else{
      norm_mat  <- matrix(
        stats::rnorm(n_bs_smp * nrow(ic_est)), nrow = n_bs_smp
        )
      e_lm_dstr <- gen_boot_sample(
        norm_mat, ic_est, center = TRUE, rate = "rootn"
        )
    }
  }
  pos_lp_norms <- control$pos_lp_norms
  nrm_type <- control$nrm_type
  num_norms <- length(pos_lp_norms)
  num_obs   <- nrow(obs_data)
  cutoff_vals <- rep(NA, num_norms)
  if (nrm_type == "bonf") {
    par_est <- sqrt(nrow(obs_data)) *
      (as.vector(apply(ic_est, 2, mean)) +
         as.vector(est_and_ic$est))
    cv_est <- NULL
    chsn_norm <- NULL
    t_s_f <- NULL
  }else{
    for (nrm_idx in 1:num_norms) {
      normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                              p = pos_lp_norms[nrm_idx], type = nrm_type)
      if (sum(is.na(normalized_obs)) != 0) browser()
      cutoff_vals[nrm_idx] <- stats::quantile(normalized_obs, 0.95)
    }
    test_stat_func <- control$test_stat_func
    # TODO: Make a helper function that does everything here until
    # the end of the if statement below: 
    if (!is.function(test_stat_func)) {
      if (test_stat_func == "est_pow") {
        t_s_f <- function(x, p) {
          pow_est <- pow_for_mag(
            boot_data = e_lm_dstr, dir = x, lp = p,
            nf_quants = cutoff_vals[which(p == pos_lp_norms)],
            nrm_type = nrm_type
            )
          return(pow_est)
          }
      } else if (test_stat_func == "mag") {
        t_s_f <- function(x, p) {
          mag_est <- mag_for_pow(
            boot_data = e_lm_dstr, dir = x,
            lp_nrms = p, power = 0.8, nrm_type = nrm_type,
            nf_quants = cutoff_vals[which(p == pos_lp_norms)]
            )
          return(mag_est)
        }
      }
    }else{
      t_s_f <- test_stat_func
    }
    ## If there is more than one row in the observed data object, it is assumed that
    ## this is the actual observed data. If there is only one row,
    ## it is assumed that this data is a single draw from the estimated
    ## limiting disitribution
    if (!exists("psi_est")) psi_est <- param_est(obs_data)$est
    par_est <- sqrt(nrow(obs_data)) * as.vector(psi_est)
    if (nrow(obs_data) > 1 & view_IC) {
       print(look_IC(ic_est))
    }
    cv_est <- as.numeric(par_est)
    test_stat <- no_fold_est(lm_dst_est = e_lm_dstr, par_est = par_est,
                             test_stat_func = t_s_f, perf_meas = perf_meas,
                             null_quants = cutoff_vals, norms_indx = pos_lp_norms,
                             norm_type = nrm_type)
    est_stats <- test_stat[[1]]
    chsn_norm <- test_stat[[2]]
    cv_est <- est_stats
  }
  if (return_lmd) {
    param_sds <- apply(
      ic_est, 2, FUN = function(x){sqrt(sum(x ** 2) / length(x))}
      )
    if (!exists("psi_est")) psi_est <- NULL
    if (exists("est_and_ic")) {
      oth_ic_inf <- est_and_ic[setdiff(names(est_and_ic),
                                       c("est", "ic"))]
    }else{
      oth_ic_inf <- NULL
    }
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm,
                "est_lm_dstr" = e_lm_dstr, "t_s_f" = t_s_f,
                "ic_est" = ic_est, "param_ests" = psi_est,
                "param_sds" = param_sds, "oth_ic" = oth_ic_inf))
  }else{
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm))
  }
}
