find_mag <- function(one_obs, dir, cutoff, nrm_idx, nrm_type){
  if(as.numeric(nrm_idx) %% 2 == 0 && nrm_type == "lp"){
    lp <- as.numeric(nrm_idx)
    n_covs <- length(one_obs)
    poly <- rep(NA, lp + 1)
    dir_mat <- matrix(NA, nrow = lp + 1, ncol = n_covs)
    dir_mat[1, ] <- rep(1, n_covs)
    for(i in 2:(lp + 1)){dir_mat[i, ] <- dir_mat[i - 1, ] * dir}
    for(k_idx in 0:lp){
      poly[k_idx + 1] <- choose(lp, k_idx) * 
        sum(dir_mat[k_idx + 1, ] * (one_obs ** (lp - k_idx)))
    }
    poly[1] <- poly[1] - cutoff ** lp
    all_rts <- polyroot(poly)
    # print(min(abs(Im(all_rts))))
    pos_root <- which(Re(all_rts) > 0)
    apr <- all_rts[pos_root]
    real_root <- which(abs(Im(apr)) < 0.00001)
    if(length(real_root) != 1){ return(-1) }else{ return(Re(apr[real_root])) }
  }else{
    obj_func <- function(c) l_p_norm(x = one_obs + c * dir, 
                                     p = nrm_idx, type = nrm_type) - cutoff
    ul <- 2
    while(obj_func(ul) <= 0){
      ul <- ul * 2
      if(abs(ul) >= 1000){stop("Cant find upper limit")}
    }
    if(obj_func(0) >= 0){
      return(-1)
    }else{
      return(uniroot(f = obj_func, lower = 0, upper = ul)$root)
    }
  }
}