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
  num_folds <- control$num_folds
  n_bs_smp <- control$n_bs_smp
  perf_meas <- control$perf_meas
  if (!is.null(control$est_lm_distr)) {
    e_lm_dstr <- est_lm_distr
  }else{
    est_and_ic <- param_est(obs_data, what = "both")
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
      if (sum(is.na(normalized_obs)) != 0) {browser()}
      cutoff_vals[nrm_idx] <- stats::quantile(normalized_obs, 0.95)
    }
    test_stat_func <- control$test_stat_func
    if (!is.function(test_stat_func)) {
      if (test_stat_func == "est_pow") {
        t_s_f <- function(x, p){
          pow_est <- pow_for_mag(
            boot_data = e_lm_dstr, dir = x, lp = p,
            nf_quants = cutoff_vals[which(p == pos_lp_norms)],
            nrm_type = nrm_type
            )
          return(pow_est)
          }
      }else if (test_stat_func == "mag") {
        t_s_f <- function(x, p){
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
    cv_est_idx <- sample(1:num_obs, replace = FALSE)
    if (num_folds == 1) {
      ## If there is more than one row in the observed data object, it is assumed that
      ## this is the actual observed data.
      ## If there is only one row, it is assumed that this data is a single draw
      ## from the estimated limiting disitribution
      if (!exists("psi_est")) {psi_est <- param_est(obs_data)$est}
      if (nrow(obs_data) > 1) {
        ## Move some of this code inside of the DE2 code.
        cor_term <- as.vector(apply(ic_est, 2, mean))
        uncor_term <- as.vector(psi_est)
        par_est <- sqrt(nrow(obs_data)) * (cor_term + uncor_term)
        cor_term <- sqrt(nrow(obs_data)) * cor_term
        if (view_IC) print(look_IC(ic_est))
      }else{
        par_est <- sqrt(nrow(obs_data)) *
          as.vector(psi_est)
      }
      if (control$nrmlize & is.null(est_lm_distr)) {
        vcov_mat <- t(ic_est) %*% ic_est / nrow(obs_data)
        # print(vcov_mat)
        if (nrow(obs_data) < nrow(vcov_mat)) {
          s_mat <- MASS::ginv(vcov_mat)
        }else{
          s_mat <- solve(vcov_mat)
        }
        # print(s_mat)
        sqrt_mat <- Re(expm::sqrtm(s_mat))
        par_est <- sqrt_mat %*% par_est
      }
      cv_est <- as.numeric(par_est)
      test_stat <- no_fold_est(lm_dst_est = e_lm_dstr, par_est = par_est,
                               test_stat_func = t_s_f, perf_meas = perf_meas,
                               null_quants = cutoff_vals,
                               norms_indx = pos_lp_norms, norm_type = nrm_type)
      est_stats <- test_stat[[1]]
      chsn_norm <- test_stat[[2]]
      cv_est <- est_stats
    }else{
      cv_est_stats <- chsn_norm <- rep(NA, num_folds)
      for (fld_idx in 1:num_folds) {
        fld_obs_idx <- c(round(num_obs * (fld_idx - 1)/num_folds) + 1,
                         round(num_obs * fld_idx/num_folds))
        training_index <- cv_est_idx[
          (fld_obs_idx[1]:fld_obs_idx[2]) * control$train_mlt
          ]
        fold_estimate <- one_fold_est(
          lm_dst_est = e_lm_dstr, obs_data = obs_data,
          param_est = param_est, test_stat_func = t_s_f,
          trn_indx = training_index, null_quants = cutoff_vals,
          norms_indx = pos_lp_norms, norm_type = nrm_type,
          perf_meas = perf_meas)
        cv_est_stats[fld_idx] <- as.numeric(fold_estimate[1])
        chsn_norm[fld_idx] <- fold_estimate[2]
      }
      cv_est <- control$f_cv_summary(cv_est_stats)
    }
  }
  if (return_lmd) {
    param_sds <- apply(
      ic_est, 2, FUN = function(x){sqrt(sum(x ** 2) / length(x))}
      )
    if (!exists("cor_term")) {cor_term <- NULL}
    if (!exists("psi_est")) {psi_est <- NULL}
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm,
                "est_lm_dstr" = e_lm_dstr, "t_s_f" = t_s_f,
                "ic_est" = ic_est, "param_ests" = psi_est,
                "param_sds" = param_sds, "correction" = cor_term))
  }else{
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm))
  }
}
