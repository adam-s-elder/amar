#' Runs a single test, with the posibility of using
#' either an apriori chosen l_p norm, or choosing
#' the norm adaptively
#'
#' @param obs_data the observed data, with columns representing covariates and
#' rows representing observations.  The first column should be the outcome of interest
#' @param lp_norm the l_p norm to be used.  Select "ap" for adaptive selection of the
#' l_p norm.
#' @param pos_lp the possible norms to be selected from.
#' @param as_meas the association measure to be used.
#' @param num_boot number of bootstrap samples used to approximate the limiting distribution
#' @param fixed the parameter (either power or magnitude) to be held constant in the
#' l_p norm selection proceedure.
#' @param cmp_meas the comparitave measure to use to select l_p norm (either mean or
#' worst direction).
#'
#'
#' @return Either the average, or minimum power for the specified
#' lp norm for all local alternatives.
#'
#' @export

cel_test <- function(obs_data, lp_norm = "ap", as_meas,
                     num_boot = 1000, pos_lp = c(2, "max"),
                     fixed = "power", cmp_meas = "mean"){
  num_obs   <- nrow(obs_data)
  est_param <- est_pearson(obs_data)
  est_cov   <- est_influence_pearson(obs_data)
  norm_mat  <- matrix(stats::rnorm(num_boot * nrow(est_cov)), nrow = num_boot)
  e_lm_dstr <- gen_boot_sample(norm_mat, est_cov, center = TRUE, rate = "rootn")
  if (lp_norm == "ap"){
    lp_v <- select_lp(e_lm_dstr, pos_lp, fxd_val = fixed, compare_measure = cmp_meas)
  }else{
    lp_v <- lp_norm
  }
  lp_norm_distr <- apply(e_lm_dstr, 1, l_p_norm, p = lp_v)
  lp_est        <- l_p_norm(est_param, p = lp_v)
  # cat(lp_norm, quantile( lp_norm_distr, 0.95), " ", sqrt(num_obs) * lp_est, "\n")
  if (stats::quantile(lp_norm_distr, 0.95) < sqrt(num_obs) * lp_est){test <- 1}else{test <- 0}
  return(list("test" = test, "norm_used" = lp_v))
}

# test_one <- run_many?]_tests(0, 2)
# apply(matrix(as.numeric(test_one), ncol = 4), 2, mean)
