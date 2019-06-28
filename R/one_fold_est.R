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
#' @param obs_data The observed data to be used for finding the optimal norm (training),
#'  and finding the test statistic (testing).  Similar to above, each row is an observation and each 
#'  column corresponds to either the outcome (first column) or a covariate. 
#' @param est_func A function that will provide the vector of parameter estimates for a subset of 
#' obs_data.
#' @param test_stat_func A function that will provide the test statistic for the 
#' given fold (using the testing data),
#'  and uses the best norm (decided on using the training data).  
#' @param trn_indx The index of observations to be used in the training set.  
#' @param null_quants The 95 percent quantiles corresponding to the limiting distribution under the 
#' various norms. The ordering of these quantiles is the same as that of the norm_indx.
#' @param norms_indx The index of the norms to be considered.  For example if we use the l_p norm, 
#' norms_indx specifies the different p's to try.  
#' @param norm_type The type of norm to be used for the test.  Generally the l_p norm 
#' @param incl_chsn_norm Boolean indicating if chosen norm index should be returned.
#' @return learned test statistic for a single fold of data
#'
#' @export

one_fold_est <- function(lm_dst_est, obs_data, est_func, test_stat_func = l_p_norm,
                         perf_meas, trn_indx, null_quants, norms_indx, 
                         norm_type){
  trn_size <- length(trn_indx)
  n_obs <- nrow(obs_data)
  if (trn_size == 0 | trn_size == n_obs){
    trn_data <- tst_data <- obs_data
  }else{
    trn_data <- obs_data[ trn_indx, , drop = FALSE]
    tst_data <- obs_data[-trn_indx, , drop = FALSE]
  }
  trn_par_est  <- sqrt(nrow(trn_data)) * as.vector(est_func(trn_data))  
  if(perf_meas == "est_pow"){
    performs <- pow_for_mag(boot_data = lm_dst_est, dir =  trn_par_est,
                            nrm_type = norm_type, lp = norms_indx,
                            nf_quant = null_quants)
    best_norm <- norms_indx[which.max(performs)]
  }else if(perf_meas == "pval"){
    performs <- pval_for_mag(boot_data = lm_dst_est, dir = trn_par_est,
                             nrm_type = norm_type, lp = norms_indx)
    best_norm <- norms_indx[which.min(performs)]
  }else if(perf_meas == "mag"){
    performs <- mag_for_pow(boot_data = lm_dst_est, dir = trn_par_est,
                            lp_nrms = norms_indx, nf_quant = null_quants,
                            nrm_type = norm_type, power = 0.8)
    #if(rnorm(1) > 3.5){print(performs)}
    best_norm <- norms_indx[which.min(performs)]
  }#else if(perf_meas == "magn"){
  #   performs <- mag_for_pow(boot_data = lm_dst_est, dir = trn_par_est,
  #                           nrm_type = )
  # }
  tst_par_est <- sqrt(nrow(tst_data)) * as.vector(est_func(tst_data))
  return(c("test_stat" = test_stat_func(tst_par_est, p = best_norm),
           "norm_choice" = best_norm))
}


