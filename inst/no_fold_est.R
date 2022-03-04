#################################################
#### Author : Adam Elder
#### Date   : November 29th 2018
#### This script is the implementation of the
#### cross validation based test that would
#### optimize the l_p norm within the procedure
#################################################

#' Helper function for cv_test.  This function estimates the cross validated
#' test statistic for a single fold of data.  Optionally this function can
#' also ouput the chosen norm index chosen using the training set.
#'
#' @param lm_dst_est The limiting distribution of the vector of parameter estimates. Each row is an
#' observation, and each column corresponds to a parameter.
#' @param par_est Estimate of the parameter to be used in deciding on a norm.
#' @param test_stat_func A function that will provide the test statistic for the
#' given fold (using the testing data),
#'  and uses the best norm (decided on using the training data).
#' @param perf_meas The preferred measure to carry out evaluation of how far away the parameter
#' estimates are away from zero.
#' @param null_quants The 95 percent quantiles corresponding to the limiting distribution under the
#' various norms. The ordering of these quantiles is the same as that of the norm_indx.
#' @param norms_indx The index of the norms to be considered.  For example if we use the l_p norm,
#' norms_indx specifies the different p's to try.
#' @param norm_type The type of norm to be used for the test.  Generally the l_p norm
#' @return learned test statistic for a single fold of data
#'
#' @export

no_fold_est <- function(lm_dst_est, par_est, test_stat_func = l_p_norm,
                         perf_meas, null_quants, norms_indx, norm_type) {
  if (perf_meas == "est_pow") {
    performs <- accept_rate(boot_data = lm_dst_est, dir = par_est,
                            nrm_type = norm_type, lp = norms_indx,
                            nf_quants = null_quants)
    best_norm <- norms_indx[which.min(performs)]
  }else if (perf_meas == "pval") {
    performs <- pval_for_mag(boot_data = lm_dst_est, dir = par_est,
                             nrm_type = norm_type, lp = norms_indx)
    best_norm <- norms_indx[which.min(performs)]
  }else if (perf_meas == "mag") {
    performs <- mag_for_pow(boot_data = lm_dst_est, dir = par_est,
                            lp_nrms = norms_indx, nf_quants = null_quants,
                            nrm_type = norm_type, power = 0.8)
    best_norm <- norms_indx[which.min(performs)]
  }
  return(list("test_stat" = test_stat_func(par_est, p = best_norm),
              "norm_choice" = best_norm))
}


