#################################################
#### Author : Adam Elder
#### Date   : November 29th 2018
#### This script is the implementation of the
#### parametric resampling based test that would
#### optimize the l_p norm within the procedure
#################################################

#' Carry out a simplified version of the Zhang and Laber test.
#' @param observed_data The observed data.
#' @param ts_sims The number of draws from the test statistic distribution
#' @param ld_sims The number of draws from the limiting distribution to
#' estimate each test statistic.
#' @return A list containing:
#'
#' - pvalue: The p-value of the test
#'
#' - test_stat: The corresponding test statistic
#'
#' - test_st_eld: The corresponding estimated limiting distribution of
#' the test statistic.
#'
#' - param_ests: The estimates of the measure of association (correlation)
#'
#' - param_ses: The corresponding standard errors of the estimates
#'
#' @export
#'
#' @examples
#'
#' ZL(data.frame(z = rnorm(100), x1 = rnorm(100), x2 = rnorm(100)), 10, 10)

ZL <- function(observed_data, ts_sims, ld_sims){
  cov_mat <- stats::cov(observed_data)
  margin_vars <- apply(observed_data[, -1], 2, stats::var)
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
  chsn_nrms <- rep(0, ncol(observed_data))
  chsn_nrms[test_stat[2]] <- 1
  return(
    list(
      "pvalue" = p_val,
      "test_stat" = test_stat[1],
      "test_st_eld" = ts_est_dstr,
      "chosen_norm" = chsn_nrms,
      "param_ests" = psi_hats,
      "param_ses" = sqrt(diag(x_cor_mat))
    )
  )
}

#' Helper function for the Zhang and Laber test.
#' @param obs_data The observed data
#' @return Returns the estimated t-statistics for the Zhang and Laber tests.

get_test_stat <- function(obs_data){
  marg_vars <- apply(obs_data, 2, stats::var)
  cov_xy    <- apply(obs_data[, -1], 2, function(x) stats::cov(obs_data[, 1], x))
  var_y <- marg_vars[1]
  var_x <- marg_vars[-1]
  betas <- cov_xy/var_x
  var_ratio <- var_y/var_x
  t_stats <- sqrt(nrow(obs_data)) * betas / sqrt(var_ratio - betas ** 2)
  return(t_stats ** 2)
}

#' A helper function to estimate the power using the generated test statistic,
#' and estimated distribution of the test statistic.
#' @param tr_lm_dstr limiting distribution of the test statistic
#' @param ts_vec Vector containing the test statistic
#' @return Return the estimated power of the test based for
#' each of the possible sum of squares norms

est_pows <- function(tr_lm_dstr, ts_vec) {
  num_opts <- length(ts_vec)
  ord_ts <- cumsum(sort(ts_vec, decreasing = TRUE))
  est_p_vals <- rep(NA, num_opts)
  for(k_idx in 1:num_opts) {
    est_p_vals[k_idx] <- mean(tr_lm_dstr[, k_idx] >= ord_ts[k_idx])
  }
  return(c(min(est_p_vals), which.min(est_p_vals)))
}

