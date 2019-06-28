###########################
##### Single Bootstrap
##### Author : Adam Elder
###########################

#' @export

single_bootstrap <- function(data, ifc, sims, norms = c("max")){
  num_norms <- length(norms)
  num_obs <- nrow(data)
  my_data <- data
  param_est_vals <- est_spearman(my_data)
  ic_ests <- est_influence_spearman(my_data)
  nu_matrix <- matrix(rnorm(sims * num_obs), nrow = num_obs, ncol = sims)
  level_1_boot <- gen_boot_sample(nu_matrix, ic_ests, center = TRUE, param_est_vals)
  l_p_res <- rep(NA, num_norms)
  for(norm_indx in 1:num_norms){
    l_p_distr <- apply(level_1_boot, 2, l_p_norm, p = norms[norm_indx])
    l_p_res[norm_indx] <- as.numeric(quantile(l_p_distr, 0.95) <
                                     l_p_norm(sqrt(num_obs) * param_est_vals,
                                              p = norms[norm_indx])) # sample(l_p_distr, 1)
  }
  return(l_p_res)
}

