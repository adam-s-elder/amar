#################################################
#### Author : Adam Elder
#### Date   : November 29th 2018
#### This script is the implementation of the
#### parametric resampling based test that would
#### optimize the l_p norm within the procedure
#################################################
#' @export

laber_test <- function(obs_data, pos_lp_norms, num_folds, n_bs_smp, nrm_type = "lp"){
  ## The cross validation procedure for our observed data
  num_obs   <- nrow(obs_data)
  num_norms <- length(pos_lp_norms)
  est_param <- est_pearson(obs_data)
  est_cov   <- est_influence_pearson(obs_data)
  norm_mat  <- matrix(rnorm(n_bs_smp * nrow(est_cov)), nrow = n_bs_smp)
  e_lm_dstr <- gen_boot_sample(norm_mat, est_cov, center = TRUE, rate = "rootn")
  cutoff_vals <- rep(NA, num_norms)
  for(nrm_idx in 1:num_norms){
    normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                            p = pos_lp_norms[nrm_idx], type = nrm_type)
    cutoff_vals[nrm_idx] <- quantile(normalized_obs, 0.95)
  }
  param_est <- est_pearson(obs_data)
  lp_perf <- rep(NA, num_norms)
  magn <- sqrt(sum(param_est ** 2)) * sqrt(length(vald_idx))
  #cat("Magn =", round(magn, 2), " ")
  for(lp_idx in 1:num_norms){
    lp_perf[lp_idx] <- pow_for_mag(e_lm_dstr, dir = param_est,
                           magn, lp = pos_lp_norms[lp_idx],
                           nf_quant = cutoff_vals[lp_idx], nrm_type = nrm_type)
  }

  test_stat <- max(lp_perf)
  ## Simulating the above proceedure for new data.
  sim_ts_mat  <- matrix(rnorm(n_bs_smp * nrow(est_cov)),
                        nrow = n_bs_smp)
  f_e_lm_dstr <- gen_boot_sample(sim_ts_mat, est_cov, center = TRUE, rate = "rootn")
  boot_sims <- rep(NA, n_bs_smp)
  for(bs in 1:n_bs_smp){
    boot_p_est<-  f_e_lm_dstr[bs, ]
    lp_perf <- rep(NA, num_norms)
    b_magn <- sqrt(sum(boot_p_est ** 2))
    #if (bs %% 10 == 0){cat("Magn =", round(magn, 2), " ")}
    for(lp_idx in 1:num_norms){
      lp_perf[lp_idx] <- pow_for_mag(e_lm_dstr, dir = boot_p_est,
                             b_magn, lp = pos_lp_norms[lp_idx],
                             nf_quant = cutoff_vals[lp_idx], nrm_type = nrm_type)
    }
    boot_sims[bs] <-  max(lp_perf)

  }
  ## Comparing the simulated null distribution to our observed value
  p_val <- mean(test_stat < boot_sims)
  return(p_val)
}


bonf_test <- function(obs_data){
  num_cov <- ncol(obs_data)
  pvals <- rep(NA, num_cov - 1)
  for(cov_idx in 2:num_cov){
    cov_lm <- lm(obs_data[, 1] ~ obs_data[, cov_idx])
    pvals[cov_idx - 1] <- summary(cov_lm)$coefficients[2, 4]
  }
  return(pvals)
}







