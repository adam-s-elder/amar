# setwd("~/Dropbox/Adam_Project/Code/L_P_Norm")
#
# source("Test_Statistic.R") #Gives the various forms of tests (one sided vs two sided)
#
# source("Est_Covariance.R") # Gives the formula that allow for calculation of the
#                            # estimates of a parameter, and the corresponding
#                            # estimated covariance matrix for the estimates
#
# source("Gen_Cop.R")        # This gives a way of generating dependent, non - normal data
#                            # The function is called gen_data
# library("MASS")
# library("methods")
# library("expm")

l_p_norm <- function(x, p = "max"){
  if (p == "max") {
    return(max(abs(x)))
  } else {
    l_p <- as.integer(p)
    return(sum(abs(x)**l_p)**(1 / l_p))
  }
}

test_one_lp_2 <- function(obs_data, norms = c("max")){
  num_obs <- nrow(obs_data)
  my_data <- obs_data
  param_est_vals <- est_spearman(my_data)
  ic_ests <- est_influence_spearman(my_data)
  mult_boot <- gen_boot_sample(1000, ic_ests, rep(0, length(param_est_vals)))
  l_p_distrs <- matrix(NA, nrow = length(norms), ncol = 1000)

  for(i in 1:length(norms)){
    l_p_distrs[i, ] <- apply(mult_boot, 1, l_p_norm, p = norms[i])
  }

  crit_95_vals <- apply(l_p_distrs, 1, function(x) stats::quantile(x, 0.95) * sqrt(num_obs))
  test_stat_vals <- vapply(X = norms, FUN = l_p_norm, FUN.VALUE = -99, x = param_est_vals, USE.NAMES = FALSE)
  test_vec <- as.numeric(test_stat_vals > crit_95_vals)
  names(test_vec) <- norms
  return(test_vec)
}

test_many_lp_2 <- function(obs_data, sims, norms = c("max")){
  total <- matrix(NA, nrow = sims, ncol = length(norms))
  ptm <- proc.time()
  est_ic <- est_influence_spearman(obs_data)
  est_param <- est_spearman(obs_data)
  for(i in 1:sims){
    boot_data <- gen_boot_sample(100, est_ic, est_param)
    total[i, ] <- test_one_lp_2(boot_data, norms = norms)
  }
  total_time <- proc.time() - ptm
  results <- c(apply(total, 2, mean,  na.rm = TRUE), total_time[3])
  names(results) <- c(norms, "Time")
  return(results)
}

find_best_p_power <- function(train_data, sims, norms = c("max")){
  find_pows <- test_many_lp_2(train_data, sims = sims, norms = norms)
  pows_no_time <- as.numeric(find_pows[ - (length(norms ) +1) ])
  return(norms[which.max( pows_no_time)  ] )
}

find_best_p_pval <- function(train_data, norms = c("max")){

  num_obs <- nrow(train_data)

  my_data <- train_data[sample(1:num_obs, replace = TRUE), ]

  param_est_vals <- est_spearman(my_data)
  ic_ests <- est_influence_spearman(my_data)

  mult_mat <- matrix(stats::rnorm(1000 * num_obs), nrow = 1000, ncol = num_obs)
  mult_boot <- mult_mat %*% ic_ests / num_obs
  l_p_distrs <- matrix(NA, nrow = length(norms), ncol = 1000)

  for(i in 1:length(norms)){
    l_p_distrs[i, ] <- apply(mult_boot, 1, l_p_norm, p = norms[i])
  }

  test_stat_vals <- vapply(X = norms, FUN = l_p_norm,
                           FUN.VALUE = -99, x = param_est_vals,
                           USE.NAMES = FALSE)

  num_norms <- length(norms)
  p_vals <- rep(NA, length = num_norms)

    for(i in 1:num_norms){
      p_vals[i] <- mean(as.numeric(l_p_distrs[i, ] <= test_stat_vals[i]))
    }

  return(norms[which.min(p_vals)])
}

choose_p_funs <- list("power" = find_best_p_power, "pval" = find_best_p_pval)

test_once_double <- function(all_data, norms = c("max"), pmethod = "power"){
  num_obs <- nrow(all_data)
  find_p_indx <- sample(1:num_obs, floor(num_obs/2), replace = FALSE)
  find_p_data <- all_data[find_p_indx, ]
  if (pmethod == "power") {
    chosen_p <- find_best_p_power(find_p_data, sims = 100, norms = norms)
  }else{
    chosen_p <- find_best_p_pval(find_p_data, norms = norms)
    }
  test <- test_one_lp_2(all_data[-find_p_indx, ], norms = chosen_p)
  return(test)
}

check_double <- function(Sigma, norms = c("max"), pmethod = "power", sims, ss){
  ptm <- proc.time()
  obs_test <- matrix(NA, ncol = length(norms), nrow = sims)
  dim <- nrow(Sigma)
  for(i in 1:sims){
    data_i <- MASS::mvrnorm(ss * 2, mu = rep(0, dim), Sigma = Sigma)
    obs_test[i, ] <- test_once_double(data_i, norms = norms, pmethod = pmethod)
  }
  total_time <- proc.time() - ptm
  pows <- c(apply(obs_test, 2, mean), total_time[3])
  names(pows) <- c(norms, "time")
  return(pows)
}

check_single <- function(Sigma, norms = c("max"), sims, ss){
  ptm <- proc.time()
  obs_test <- matrix(NA, ncol = length(norms), nrow = sims)
  dim <- nrow(Sigma)
  for(i in 1:sims){
    data_i <- MASS::mvrnorm(ss, mu = rep(0, dim), Sigma = Sigma)
    obs_test[i, ] <- test_one_lp_2(data_i, norms = norms)
  }
  total_time <- proc.time() - ptm
  pows <- c(apply(obs_test, 2, mean), total_time[3])
  names(pows) <- c(norms, "time")
  return(pows)
}
