#' CV-test runs a cross-validation, parametric bootstrap test. This function returns
#' an approximate p-value for the specified test statistic.
#'
#' @param obs_data The observed data to be used for finding the optimal norm (training),
#' and finding the test statistic (testing).  Similar to above, each row is an observation and each
#' column corresponds to either the outcome (first column) or a covariate.
#' @param pos_lp_norms The index of the norms to be considered.  For example if we use the l_p norm,
#' norms_indx specifies the different p's to try.
#' @param num_folds The number of folds to be used in the cross-validation proceedure.  If set to 1,
#' no cross validation will be used.
#' @param f_cv_summary How test statistics from different folds are combined to create an overall
#' test statistic.  Usually the mean is used.
#' @param n_bs_smp Number of samples to be used in estimating the limiting distribution of the
#' test statistic under the null.
#' @param nrm_type The type of norm to be used for the test.  Generally the l_p norm
#' @param big_train Data is split into splits of roughly equal sizes. The number of splits is equal
#' to num_folds. If big_train is TRUE then all but one of these splits will be training data,
#' if big_train is FALSE all but one will be testing data.
#' @param test_stat_func A function that will provide the test statistic for the
#' given fold (using the testing data), and uses the best
#' norm (decided on using the training data).
#' @param incl_chsn_norm Boolean indicating if chosen norm index should be returned.
#' @return learned test statistic for a single fold of data
#'
#' @export

cv_test <- function(obs_data, param, pos_lp_norms, num_folds, f_cv_summary = mean,
                    n_bs_smp, nrm_type = "lp", test_stat_func,
                    ts_ld_bs_samp = 250, big_train = TRUE, incl_chsn_norm = FALSE,
                    show_hist = FALSE, test_type = "par_boot", perf_meas = "est_pow"){
  train_mlt <- (-1) ** (2 + as.integer(big_train))
  est_infl_func <- get_infl[[param]]$est_IC
  est_parm_fun <- f_estimate <- get_infl[[param]]$est_param
  init_cv_est <- est_cv(obs_data = obs_data, est_infl = est_infl_func,
                        pos_lp_norms = pos_lp_norms, nrm_type = nrm_type, est_func = f_estimate,
                        n_bs_smp = n_bs_smp, test_stat_func = test_stat_func,
                        f_cv_summary = f_cv_summary, trn_mlt = train_mlt,
                        num_folds = num_folds, pref_meas = perf_meas,
                        est_lm_distr =  NULL, return_lmd = TRUE)
  e_lm_dstr <- init_cv_est$est_lm_dstr
  t_s_f <- init_cv_est$t_s_f
  cv_est <- init_cv_est$cv_est
  est_cov <- init_cv_est$est_co
  if (test_type == "par_boot"){
    sim_ts_mat  <- matrix(rnorm(num_folds * ts_ld_bs_samp * nrow(est_cov)),
                          nrow = num_folds * ts_ld_bs_samp)
    f_e_lm_dstr <- gen_boot_sample(sim_ts_mat, est_cov, center = TRUE, rate = "rootn")
    ts_lim_dist <- rep(NA, ts_ld_bs_samp)
    for(bs_idx in 1:ts_ld_bs_samp){
      sub_data <- f_e_lm_dstr[(num_folds * (bs_idx - 1) + 1):(num_folds * bs_idx), ,drop = FALSE]
      par_boot_cv_est <- est_cv(obs_data = sub_data, est_infl = est_infl_func,
                                pos_lp_norms = pos_lp_norms, nrm_type = nrm_type,
                                n_bs_smp = n_bs_smp, test_stat_func = t_s_f,
                                f_cv_summary = f_cv_summary,
                                est_func = function(x)apply(x, 2, mean), pref_meas = perf_meas,
                                trn_mlt = train_mlt, num_folds = num_folds,
                                est_lm_distr = e_lm_dstr, return_lmd = FALSE)
      ts_lim_dist[bs_idx] <- par_boot_cv_est$cv_est
    }
  }else if(test_type == "perm"){
    ts_lim_dist <- rep(NA, ts_ld_bs_samp)
    num_obs <- nrow(obs_data)
    for(perm_idx in 1:ts_ld_bs_samp){
      y_idx <- sample(1:num_obs, replace = FALSE)
      perm_data <- obs_data
      perm_data[, 1] <- perm_data[y_idx, 1]
      perm_cv_est <- est_cv(obs_data = perm_data, est_infl = est_infl_func,
                            pos_lp_norms = pos_lp_norms, nrm_type = nrm_type,
                            f_cv_summary = f_cv_summary,
                            est_func = f_estimate, n_bs_smp = n_bs_smp, pref_meas = perf_meas,
                            test_stat_func = test_stat_func, trn_mlt = train_mlt,
                            num_folds = num_folds, est_lm_distr = NULL, return_lmd = FALSE)
      ts_lim_dist[perm_idx] <- perm_cv_est$cv_est
    }
  }
  if (show_hist){hist(ts_lim_dist);abline(v = cv_est, col = "red");browser()}
  init_p_val <- mean(as.integer(cv_est <= ts_lim_dist))
  if(test_stat_func %in% c("mag", "pval")){
    p_val <- 1 - init_p_val
  }else{
    p_val <- init_p_val
  }
  if (incl_chsn_norm){
    chsn_tbl <- vapply(pos_lp_norms, function(x) mean(x == init_cv_est$chsn_norms),
                       FUN.VALUE = -99)
    return(c("Estiamted p-value" = p_val,
             chsn_tbl))
  }else{
    return(c("Estimate p-value" = p_val))
  }
}

