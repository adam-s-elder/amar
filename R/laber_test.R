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

#' @export
ZL <- function(observed_data, ts_sims, ld_sims){
  cov_mat <- cov(observed_data)
  margin_vars <- apply(observed_data[, -1], 2, var)
  x_cor_mat <- diag(1/sqrt(margin_vars)) %*%
    cov_mat[-1, -1] %*%
    diag(1/sqrt(margin_vars))
  psi_hats <- get_test_stat(observed_data)
  null_lm_distr <- MASS::mvrnorm(n = ts_sims, mu = rep(0, length(psi_hats)),
                                 Sigma = x_cor_mat)
  null_trs_dstr <- matrix(NA, nrow = nrow(null_lm_distr), ncol = ncol(null_lm_distr))
  for(row_idx in 1:nrow(null_lm_distr)){
    sim_row <- null_lm_distr[row_idx, ] ** 2
    null_trs_dstr[row_idx, ] <- cumsum(sort(sim_row, decreasing = TRUE))
  }
  test_stat <- est_pows(null_trs_dstr, psi_hats)
  draws <- MASS::mvrnorm(n = ld_sims, mu = rep(0, length(psi_hats)),
                         Sigma = x_cor_mat)
  ts_est_dstr <- rep(NA, ld_sims)
  for(d_idx in 1:ld_sims){
    draw <- draws[d_idx, ]
    d_ts <- draw ** 2
    ts_est_dstr[d_idx] <- est_pows(null_trs_dstr, d_ts)[1]
  }
  p_val <- mean(test_stat[1] >= ts_est_dstr)
  return(p_val)
}

#' @export
get_test_stat <- function(obs_data){
  marg_vars <- apply(obs_data, 2, var)
  cov_xy    <- apply(obs_data[, -1], 2, function(x) cov(obs_data[, 1], x))
  var_y <- marg_vars[1]
  var_x <- marg_vars[-1]
  betas <- cov_xy/var_x
  var_ratio <- var_y/var_x
  t_stats <- sqrt(nrow(obs_data)) * betas / sqrt(var_ratio - betas ** 2)
  return(t_stats ** 2)
}

#' @export
est_pows <- function(tr_lm_dstr, ts_vec){
  num_opts <- length(ts_vec)
  ord_ts <- cumsum(sort(ts_vec, decreasing = TRUE))
  est_p_vals <- rep(NA, num_opts)
  for(k_idx in 1:num_opts){
    est_p_vals[k_idx] <- mean(tr_lm_dstr[, k_idx] >= ord_ts[k_idx])
  }
  return(c(min(est_p_vals), which.min(est_p_vals)))
}

#' @export
bonf_test <- function(obs_data, test_type = "pearson"){
  if(test_type == "pearson"){
    num_cov <- ncol(obs_data)
    pvals <- rep(NA, num_cov - 1)
    for(cov_idx in 2:num_cov){
      cov_lm <- lm(obs_data[, 1] ~ obs_data[, cov_idx])
      pvals[cov_idx - 1] <- summary(cov_lm)$coefficients[2, 4]
    }
  }ifelse(test_type == "mean"){
    num_cov <- ncol(obs_data)
    pvals <- rep(NA, num_cov)
    for(cov_idx in 1:num_cov){
      pvals[cov_idx] <- t.test(obs_data[, cov_idx])$p.value
    }
  }
  return(pvals)
}


#' @export
ZL_use_infl <- function(observed_data, ts_sims, ld_sims){
  cov_mat <- cov(observed_data)
  margin_vars <- apply(observed_data[, -1], 2, var)
  x_cor_mat_p <- est_influence_pearson(observed_data)
  x_cor_mat <- t(x_cor_mat_p) %*% x_cor_mat_p / nrow(observed_data)
  psi_hats <- get_test_stat(observed_data)
  null_lm_distr <- MASS::mvrnorm(n = ts_sims, mu = rep(0, length(psi_hats)),
                                 Sigma = x_cor_mat)
  null_trs_dstr <- matrix(NA, nrow = nrow(null_lm_distr),
                          ncol = ncol(null_lm_distr))
  for(row_idx in 1:nrow(null_lm_distr)){
    sim_row <- null_lm_distr[row_idx, ] ** 2
    null_trs_dstr[row_idx, ] <- cumsum(sort(sim_row, decreasing = TRUE))
  }
  test_stat <- est_pows(null_trs_dstr, psi_hats)
  draws <- MASS::mvrnorm(n = ld_sims, mu = rep(0, length(psi_hats)),
                         Sigma = x_cor_mat)
  ts_est_dstr <- rep(NA, ld_sims)
  for(d_idx in 1:ld_sims){
    draw <- draws[d_idx, ]
    d_ts <- draw ** 2
    ts_est_dstr[d_idx] <- est_pows(null_trs_dstr, d_ts)[1]
  }
  p_val <- mean(test_stat[1] >= ts_est_dstr)
  return(p_val)
}




