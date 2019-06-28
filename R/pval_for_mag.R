pval_for_mag <- function(boot_data, dir, lp = 2, nrm_type = "lp"){
  num_norms <- length(lp)
  est_pvals <- rep(NA, num_norms)
  
  distr <- matrix(NA, nrow = nrow(boot_data), ncol = num_norms)
  if(nrm_type == "ordl2"){
    num_obs <- nrow(boot_data)
    trans_dir <- cumsum(dir ** 2, decreasing = TRUE)
    for(obs_idx in 1:num_obs){
      distr[obs_idx, ] <- cumsum(sort(boot_data[obs_idx, ] ** 2, decreasing = TRUE))
    }
    for(lp_idx in 1:num_norms){
      trans_est <- trans_dir[lp[lp_idx]]
      est_pvals[lp_idx] <- mean(as.numeric(distr[, lp[lp_idx]] >= trans_est))
    }
    return(est_pvals)
  }
  if(nrm_type == "lp"){
    nrmd_distr <- apply(boot_data, 1, l_p_norm, p = lp[lp_idx], type = "lp")
    for(lp_idx in 1:num_norms){
      trans_dir <- l_p_norm(dir, p = lp[lp_idx], type = "lp")
      est_pvals[lp_idx] <- mean(as.numeric(nrmd_distr >= trans_dir))
    }
    return(est_pvals)
  }
}