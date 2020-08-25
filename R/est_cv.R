#' This is an intermediate function inbetween the cv_test and one_fold_est functions which
#' allows for specification of the measure of performance.
#'
#' @param obs_data The observed data to be used for finding the optimal norm (training),
#' and finding the test statistic (testing).  Similar to above, each row is an observation and each
#' column corresponds to either the outcome (first column) or a covariate.
#' @param est_infl Estimate of the influence function, if one exists.
#' @param est_func Function that estimates the parameter.
#' @param pos_lp_norms The index of the norms to be considered.  For example if we use the l_p norm,
#' norms_indx specifies the different p's to try.
#' @param num_folds The number of folds to be used in the cross-validation proceedure.  If set to 1,
#' no cross validation will be used.
#' @param f_cv_summary How test statistics from different folds are combined to create an overall
#' test statistic.  Usually the mean is used.
#' @param trn_mlt Multiplier which specifies if the larger set of data will be in the training or
#' testing set.
#' @param return_lmd Boolean for whether or not to return the estimated limiting distribution.
#' @param est_lm_distr An estimate of the limiting distribution if it is provided
#' @param pref_meas the prefered measure used to generate the test statistic.
#' @param n_bs_smp Number of samples to be used in estimating the limiting distribution of the
#' test statistic under the null.
#' @param nrm_type The type of norm to be used for the test.  Generally the l_p norm.
#' @param test_stat_func A function that will provide the test statistic for the
#' given fold (using the testing data), and uses the best
#' norm (decided on using the training data).
#' @param nrmlize Boolean for if the estimator should be normalized to have a limiting distribution with identity covaraince matrix.
#' @return learned test statistic for a single fold of data
#'
#' @export

est_cv <- function(obs_data, est_infl, pos_lp_norms, nrm_type, est_func,
                   n_bs_smp, test_stat_func, trn_mlt, num_folds, pref_meas,
                   f_cv_summary = mean, est_lm_distr = NULL,
                   return_lmd = FALSE, view_IC = FALSE,
                   nrmlize = FALSE){
  if (!is.null(est_lm_distr)){
    e_lm_dstr <- est_lm_distr
  }else{
    est_cov <- est_infl(obs_data)
    if(nrmlize){
      e_lm_dstr <- matrix(stats::rnorm(n_bs_smp * ncol(est_cov)), nrow = n_bs_smp)
    }else{
      norm_mat  <- matrix(stats::rnorm(n_bs_smp * nrow(est_cov)), nrow = n_bs_smp)
      e_lm_dstr <- gen_boot_sample(norm_mat, est_cov,
                                   center = TRUE, rate = "rootn")
    }
  }
  num_norms <- length(pos_lp_norms)
  num_obs   <- nrow(obs_data)
  cutoff_vals <- rep(NA, num_norms)
  if(nrm_type == "bonf"){
    param_est <- sqrt(nrow(obs_data)) *
      (as.vector(apply(est_cov, 2, mean)) +
         as.vector(est_func(obs_data)))
    cv_est <- NULL
    chsn_norm <- NULL
    t_s_f <- NULL
  }else{
    for(nrm_idx in 1:num_norms){
      normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                              p = pos_lp_norms[nrm_idx], type = nrm_type)
      if(sum(is.na(normalized_obs)) != 0){browser()}
      cutoff_vals[nrm_idx] <- stats::quantile(normalized_obs, 0.95)
    }
    if (!is.function(test_stat_func)){
      if (test_stat_func == "est_pow"){
        t_s_f <- function(x, p){
          pow_est <- pow_for_mag(boot_data = e_lm_dstr, dir = x, lp = p,
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
    if(num_folds == 1){
      ## If there is more than one row in the observed data object, it is assumed that
      ## this is the actual observed data.
      ## If there is only one row, it is assumed that this data is a single draw
      ## from the estimated limiting disitribution
      if (nrow(obs_data) > 1) {
        cor_term <- as.vector(apply(est_cov, 2, mean))
        uncor_term <- as.vector(est_func(obs_data))
        param_est <- sqrt(nrow(obs_data)) * (cor_term + uncor_term)
        cor_term <- sqrt(nrow(obs_data)) * cor_term
        if (view_IC) print(look_IC(est_cov))
      }else{
        param_est <- sqrt(nrow(obs_data)) * as.vector(est_func(obs_data))
      }
      if(nrmlize & is.null(est_lm_distr)){
        vcov_mat <- t(est_cov) %*% est_cov / nrow(obs_data)
        # print(vcov_mat)
        if(nrow(obs_data) < nrow(vcov_mat)){
          s_mat <- MASS::ginv(vcov_mat)
        }else{
          s_mat <- solve(vcov_mat)
        }
        # print(s_mat)
        sqrt_mat <- Re(expm::sqrtm(s_mat))
        param_est <- sqrt_mat %*% param_est
      }
      cv_est <- as.numeric(param_est)
      test_stat <- no_fold_est(lm_dst_est = e_lm_dstr, par_est = param_est,
                               test_stat_func = t_s_f, perf_meas = pref_meas,
                               null_quants = cutoff_vals, norms_indx = pos_lp_norms,
                               norm_type = nrm_type)
      est_stats <- test_stat[[1]]
      chsn_norm <- test_stat[[2]]
      cv_est <- est_stats
    }else{

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
    }
  }
  if (return_lmd) {
    param_sds <- apply(est_cov, 2,
                      FUN = function(x){
                        sqrt(sum(x ** 2) / length(x))
                      })
    if (!exists("cor_term")) {cor_term <- NULL}
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm,
                "est_lm_dstr" = e_lm_dstr, "t_s_f" = t_s_f,
                "est_cov" = est_cov, "param_ests" = param_est,
                "param_sds" = param_sds, "correction" = cor_term))
  }else{
    return(list("cv_est" = cv_est, "chsn_norms" = chsn_norm))
  }
}
