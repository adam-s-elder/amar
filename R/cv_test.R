#################################################
#### Author : Adam Elder
#### Date   : November 29th 2018
#### This script is the implementation of the
#### cross validation based test that would
#### optimize the l_p norm within the procedure
#################################################
#################################################
#' @export



one_fold_est <- function(lm_dst_est, obs_data, est_func, sum_func,
                         trn_indx, null_quants, norms_indx, norm_type){
  # num_obs <- nrow(obs_data)
  trn_data <- obs_data[-trn_indx,]
  tst_data <- obs_data[ trn_indx,]
  #cat("train_Rows =", nrow(trn_data), ", ",
  #   "test_rows =", nrow(tst_data), "\n")
  if(nrow(trn_data) == 0){trn_data <- obs_data}
  trn_par_est  <- est_func(trn_data)
  performs <- pow_for_mag(boot_data = lm_dst_est, dir = trn_par_est,
                          nrm_type = norm_type, lp = norms_indx,
                          nf_quant = null_quants)
  best_norm <- norms_indx[which.max(performs)]
  if (nrow(tst_data) > 0){
    tst_par_est <- est_func(tst_data)
  }else{tst_par_est <- trn_par_est}
  return(sum_func(tst_par_est, p = best_norm))
}

one_fold_sim <- function(lm_dst_est, sim_ests, sum_func, norm_type,
                         trn_indx, null_quants, norms_indx){
  num_idx <- length(trn_indx)
  num_sim_ests <- nrow(sim_ests)
  if(length(trn_indx) == 0 || length(trn_indx) == num_sim_ests){
    trn_indx <- -1 * 1:num_sim_ests
    num_idx <- length(trn_indx)}
  trn_par_est  <- apply(sim_ests[-trn_indx, , drop = FALSE], 2, mean) * sqrt(num_idx)
  performs <- pow_for_mag(boot_data = lm_dst_est, dir = trn_par_est,
                          nrm_type = norm_type, lp = norms_indx,
                          nf_quant = null_quants)
  best_norm <- norms_indx[which.max(performs)]
  if (num_idx != num_sim_ests){
    tst_par_est <- apply(sim_ests[ trn_indx, , drop = FALSE], 2, mean) *
      sqrt(nrow(sim_ests) - num_idx)
  }else{tst_par_est <- trn_par_est}
  return(sum_func(tst_par_est, p = best_norm))
}


cv_test <- function(obs_data, pos_lp_norms, num_folds, f_cv_summary = mean,
                    n_bs_smp, nrm_type = "lp", big_train = TRUE, f_summary,
                    f_estimate){
  train_mlt <- (-1) ** (1 + as.integer(big_train))
  num_obs   <- nrow(obs_data)
  num_norms <- length(pos_lp_norms)
  est_cov   <- est_influence_pearson(obs_data)
  norm_mat  <- matrix(rnorm(n_bs_smp * nrow(est_cov)), nrow = n_bs_smp)
  e_lm_dstr <- gen_boot_sample(norm_mat, est_cov, center = TRUE, rate = "rootn")
  cutoff_vals <- rep(NA, num_norms)
  for(nrm_idx in 1:num_norms){
    normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                            p = pos_lp_norms[nrm_idx], type = nrm_type)
    cutoff_vals[nrm_idx] <- quantile(normalized_obs, 0.95)
  }
  cv_est_idx <- sample(1:num_obs, replace = FALSE)
  cv_est_stats <- rep(NA, num_folds)
  for(fld_idx in 1:num_folds){
    fld_obs_idx <- c(round(num_obs * (fld_idx - 1)/num_folds) + 1,
                     round(num_obs * fld_idx/num_folds))
    training_index <- cv_est_idx[(fld_obs_idx[1] : fld_obs_idx[2]) * train_mlt]
    cv_est_stats[fld_idx] <- one_fold_est(lm_dst_est = e_lm_dstr, obs_data = obs_data,
                                          est_func = f_estimate, sum_func = f_summary,
                                          trn_indx = training_index,
                                          null_quants = cutoff_vals,
                                          norms_indx = pos_lp_norms, norm_type = nrm_type)
  }
  cv_est <- f_cv_summary(cv_est_stats)
  sim_ts_mat  <- matrix(rnorm(num_folds * n_bs_smp * nrow(est_cov))/sqrt(num_obs/num_folds),
                        nrow = num_folds * n_bs_smp)
  f_e_lm_dstr <- gen_boot_sample(sim_ts_mat, est_cov, center = TRUE, rate = "rootn")
  mc_draw <- rep(NA, n_bs_smp)
  for(bs_idx in 1:n_bs_smp){
    sub_data <- f_e_lm_dstr[(num_folds * (bs_idx - 1) + 1):(num_folds * bs_idx), ,drop = FALSE]
    sim_cv_est <- rep(NA, num_folds)
    for(fld_idx in 1:num_folds){
      training_index <- fld_idx * train_mlt
      sim_cv_est[fld_idx] <- one_fold_sim(lm_dst_est = e_lm_dstr, sim_ests = sub_data,
                                          sum_func = f_summary, norm_type = nrm_type,
                                          trn_indx = training_index, null_quants = cutoff_vals,
                                          norms_indx = pos_lp_norms)
    }
    mc_draw[bs_idx] <- f_cv_summary(sim_cv_est)
  }
  return(mean(as.integer(cv_est <= mc_draw)))
}

perm_test <- function(obs_data, pos_lp_norms, num_folds, f_cv_summary = mean,
                      num_perms, nrm_type = "lp", big_train = TRUE, f_summary,
                      f_estimate, n_bs_smp = 1000, seed = NULL){
  train_mlt <- (-1) ** (1 + as.integer(big_train))
  num_obs   <- nrow(obs_data)
  num_norms <- length(pos_lp_norms)
  if(is.null(seed)){set.seed(sample(1:1000, 1))}else{set.seed(seed)}
  ## Initial Estimate
  est_cov   <- est_influence_pearson(obs_data)
  norm_mat  <- matrix(rnorm(n_bs_smp * nrow(est_cov)), nrow = n_bs_smp)
  e_lm_dstr <- gen_boot_sample(norm_mat, est_cov, center = TRUE, rate = "rootn")
  cutoff_vals <- rep(NA, num_norms)
  for(nrm_idx in 1:num_norms){
    normalized_obs <- apply(e_lm_dstr, 1, l_p_norm,
                            p = pos_lp_norms[nrm_idx], type = nrm_type)
    cutoff_vals[nrm_idx] <- quantile(normalized_obs, 0.95)
  }
  cv_est_idx <- sample(1:num_obs, replace = FALSE)
  cv_est_stats <- rep(NA, num_folds)
  for(fld_idx in 1:num_folds){
    fld_obs_idx <- c(round(num_obs * (fld_idx - 1)/num_folds) + 1,
                     round(num_obs * fld_idx/num_folds))
    training_index <- cv_est_idx[(fld_obs_idx[1] : fld_obs_idx[2]) * train_mlt]
    cv_est_stats[fld_idx] <- one_fold_est(lm_dst_est = e_lm_dstr, obs_data = obs_data,
                                          est_func = f_estimate, sum_func = f_summary,
                                          trn_indx = training_index,
                                          null_quants = cutoff_vals,
                                          norms_indx = pos_lp_norms, norm_type = nrm_type)
  }
  cv_est <- f_cv_summary(cv_est_stats)

  ## Running the permutation tesy
  perm_results <- rep(NA, num_perms)
  for(perm_idx in 1:num_perms){
    y_idx <- sample(1:num_obs, replace = FALSE)
    perm_data <- obs_data
    perm_data[, 1] <- perm_data[y_idx, 1]
    perm_est_idx <- sample(1:num_obs, replace = FALSE)
    est_perm_cov   <- est_influence_pearson(perm_data)
    norm_perm_mat  <- matrix(rnorm(n_bs_smp * nrow(est_perm_cov)), nrow = n_bs_smp)
    e_perm_lm_dstr <- gen_boot_sample(norm_perm_mat, est_perm_cov, center = TRUE, rate = "rootn")
    cutoff_perm_vals <- rep(NA, num_norms)
    for(nrm_idx in 1:num_norms){
      normalized_perm_obs <- apply(e_perm_lm_dstr, 1, l_p_norm,
                              p = pos_lp_norms[nrm_idx], type = nrm_type)
      cutoff_perm_vals[nrm_idx] <- quantile(normalized_perm_obs, 0.95)
    }
    perm_est_stats <- rep(NA, num_folds)
    for(fld_idx in 1:num_folds){
      fld_obs_idx <- c(round(num_obs * (fld_idx - 1)/num_folds) + 1,
                       round(num_obs * fld_idx/num_folds))
      training_index <- cv_est_idx[(fld_obs_idx[1] : fld_obs_idx[2]) * train_mlt]
      perm_est_stats[fld_idx] <- one_fold_est(lm_dst_est = e_perm_lm_dstr, obs_data = perm_data,
                                              est_func = f_estimate, sum_func = f_summary,
                                              trn_indx = training_index,
                                              null_quants = cutoff_perm_vals,
                                              norms_indx = pos_lp_norms, norm_type = nrm_type)
    }
    perm_results[perm_idx] <- f_cv_summary(perm_est_stats)

  }
  return(mean(perm_results >= cv_est))
}

l_p_norm <- function(x, p = "max", type = "lp"){
  if(type == "lp"){
    if (p == "max") {
      return(max(abs(x)))
    } else {
      l_p <- as.integer(p)
      return(sum(abs(x)**l_p)**(1 / l_p))
    }
  }else if (type == "ordl2"){
    l_p <- as.integer(p)
    some_x <- x[order(x)]
    return(sum(some_x[1:p] **2))
  }
}



est_influence_pearson <- function(observ){
  n <- nrow(observ)
  num_cov <- ncol(observ) - 1
  means <- colMeans(observ)
  y_mean <- means[1]
  cent_obs <- observ - matrix(rep(colMeans(observ), each = n), nrow = n)
  cent_obs_sqrd <- cent_obs^2
  sigmas <- colSums(cent_obs_sqrd)/(n - 1)
  covs <- as.numeric(crossprod(cent_obs[, 1], cent_obs[, -1]))/(n - 1)
  y_var <- sigmas[1]
  cent_cov <- cent_obs[, -1] * cent_obs[, 1] - rep(covs, each = n)
  cent_var <- cent_obs_sqrd - rep(sigmas, each = n)
  psi_1 <- matrix(rep(sigmas[-1], each = n), nrow = n) *
    cent_var[, 1] + y_var * cent_var[, -1]
  ic <- (cent_cov - psi_1 * rep(covs/(y_var * sigmas[-1]), each = n)) *
    rep((1/sqrt(y_var * sigmas[-1])), each = n)
  return(ic)
}

est_pearson <- function(observ){
  num_var <- ncol(observ)
  return(cor(observ[, 1], observ[, -1], method = "pearson"))
}




# for(bb in   1:300){
#   if(bb %% 10 == 0){cat(bb, " ")}
#   p_val <- cv_test(mvrnorm(100, mu = rep(0, 5), diag(5)),
#                    pos_lp_norms = 1:4,
#                    num_folds = 3, n_bs_smp = 100, nrm_type = "ordl2")
#   test_res_two[bb] <- p_val
# }






