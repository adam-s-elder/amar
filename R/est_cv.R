#' This is an intermediate function inbetween the cv_test and one_fold_est functions which
#' allows for specification of the measure of performance.
#'
#' @param obs_data The observed data to be used for finding the optimal norm (training),
#' and finding the test statistic (testing).  Similar to above, each row is an observation and each
#' column corresponds to either the outcome (first column) or a covariate.
#' @param est_infl Estimate of the influence function, if one exists.
#' @param pos_lp_norms The index of the norms to be considered.  For example if we use the l_p norm,
#' norms_indx specifies the different p's to try.
#' @param num_folds The number of folds to be used in the cross-validation proceedure.  If set to 1,
#' no cross validation will be used.
#' @param f_cv_summary How test statistics from different folds are combined to create an overall
#' test statistic.  Usually the mean is used.
#' @param n_bs_smp Number of samples to be used in estimating the limiting distribution of the
#' test statistic under the null.
#' @param nrm_type The type of norm to be used for the test.  Generally the l_p norm.
#' @param big_train Data is split into splits of roughly equal sizes. The number of splits is equal
#' to num_folds. If big_train is TRUE then all but one of these splits will be training data,
#' if big_train is FALSE all but one will be testing data.
#' @param num_folds The number of folds desired for the cross-validated parameter estimate (usually 1).
#' @param test_stat_func A function that will provide the test statistic for the
#' given fold (using the testing data), and uses the best
#' norm (decided on using the training data).
#' @param incl_chsn_norm Boolean indicating if chosen norm index should be returned.
#' @return learned test statistic for a single fold of data
#'
#' @export

est_cv <- function(obs_data, est_infl, pos_lp_norms, nrm_type, est_func,
                   n_bs_smp, test_stat_func, trn_mlt, num_folds, pref_meas,
                   f_cv_summary = mean, est_lm_distr = NULL, return_lmd = FALSE){
  if (!is.null(est_lm_distr)){
    e_lm_dstr <- est_lm_distr
  }else{
    est_cov   <- est_infl(obs_data)
    norm_mat  <- matrix(rnorm(n_bs_smp * nrow(est_cov)), nrow = n_bs_smp)
    e_lm_dstr <- gen_boot_sample(norm_mat, est_cov, center = TRUE, rate = "rootn")
  }
  num_norms <- length(pos_lp_norms)
  num_obs   <- nrow(obs_data)
  cutoff_vals <- rep(NA, num_norms)
  for(nrm_idx in 1:num_norms){
    normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                            p = pos_lp_norms[nrm_idx], type = nrm_type)
    if(sum(is.na(normalized_obs)) != 0){browser()}
    cutoff_vals[nrm_idx] <- quantile(normalized_obs, 0.95)
  }
  if (!is.function(test_stat_func)){
    if (test_stat_func == "est_pow"){
      t_s_f <- function(x, p){
        pow_est <- pow_for_mag(boot_data = e_lm_dstr,dir = x, lp = p,
                               nf_quants = cutoff_vals[which(p == pos_lp_norms)],
                               nrm_type = nrm_type)
        return(pow_est) }
    }else if (test_stat_func == "mag"){
      t_s_f <- function(x, p){
        mag_est <- mag_for_pow(boot_data = e_lm_dstr,dir = x, lp_nrms = p,
                               power = 0.8, nrm_type = nrm_type,
                               nf_quants = cutoff_vals[which(p == pos_lp_norms)])
        return(mag_est)
      }
    }
  }else{ t_s_f <- test_stat_func }

  cv_est_idx <- sample(1:num_obs, replace = FALSE)
  cv_est_stats <- chsn_norm <- rep(NA, num_folds)

  for(fld_idx in 1:num_folds){
    fld_obs_idx <- c(round(num_obs * (fld_idx - 1)/num_folds) + 1,
                     round(num_obs * fld_idx/num_folds))
    training_index <- cv_est_idx[(fld_obs_idx[1] : fld_obs_idx[2]) * trn_mlt]
    fold_estimate <- one_fold_est(lm_dst_est = e_lm_dstr, obs_data = obs_data,
                                  est_func = est_func,  test_stat_func = t_s_f,
                                  trn_indx = training_index, null_quants = cutoff_vals,
                                  norms_indx = pos_lp_norms, norm_type = nrm_type,
                                  perf_meas = pref_meas)
    cv_est_stats[fld_idx] <- as.numeric(fold_estimate[1])
    chsn_norm[fld_idx] <- fold_estimate[2]
  }
  cv_est <- f_cv_summary(cv_est_stats)
  if (return_lmd){
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm,
                "est_lm_dstr" = e_lm_dstr, "t_s_f" = t_s_f,
                "est_cov" = est_cov))
  }else{
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm))
  }
}
