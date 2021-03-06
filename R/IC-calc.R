######################################################################
#### Estimating Influence Curve and parameter for various parameters.
#### Adam Elder
#### This script will contain all of the various proceedures
#### used for estimating both the parameter, and the influence
#### curves used for estimating spearman correlation,
#### pearson correlation and mean.
#### These functions will be given an entire dataset, and will
#### return the estimate for the parameter, and for each
#### estimated influence curve at each observation.
######################################################################

#' The following function will take a set of observations, and return the
# estimated covariance matrix for the estimates of the pearson correlation.
# Estimates of the covariance are generated using the empirical influence
# function.  The first column of your data should correspond to the
#' variable of interest (the variable for which pearson correlation is
#' calculated).
#' @param observ the observed data.  The first column should be the outcome.
#' @param what the desired return value. Should be one of `"ic"`
#' (infludence curve), `"est"` (estimate), or `"both"`.
#' @param control any other control parameters to be passed to the estimator.
#'
#' @export

ic.pearson <- function(observ, what = "both", control = NULL){
  if (!(what %in% c("ic", "est", "both"))) {
    stop("what must be one of ic (influence curve), est (estimate), or both")
  }
  ret <- list()
  if (what %in% c("ic", "both")) {
    n <- nrow(observ)
    num_cov <- ncol(observ) - 1
    means <- colMeans(observ)
    obs_cent <- observ - matrix(
      rep(means, each = n),
      nrow = n, byrow = FALSE
    )
    sigmas <- apply(
      obs_cent, 2,
      function(x) mean(x ** 2) ** (1 / 2)
    )
    covs <- as.numeric(crossprod(obs_cent[, 1, drop = FALSE],
                                 obs_cent[, -1]))/n
    x_sd <- sigmas[-1]
    y_sd <- sigmas[1]
    cors <- covs / (y_sd * x_sd)
    cent_cov <-
      sweep(obs_cent[, -1, drop = FALSE], MARGIN = 1, obs_cent[, 1], `*`)

    piece_one <-  sweep(cent_cov, MARGIN = 2,
                        1 / (y_sd * x_sd), `*`)
    piece_two <- sweep(obs_cent[, -1, drop = FALSE] ** 2, MARGIN = 2,
                       1 / (x_sd ** 2), `*`)
    piece_three <- matrix(
      rep(obs_cent[, 1] ** 2 / y_sd ** 2, each = num_cov),
      ncol = num_cov, byrow = TRUE
    )

    ic <- piece_one -
      sweep(piece_two + piece_three,
            MARGIN = 2, (1 / 2) * cors, `*`)

    ret$ic <- ic
  }
  if (what %in% c("both", "est") ) {
    ret$est <- stats::cor(observ[, 1], observ[, -1], method = "pearson")
  }
  return(ret)
}

## Check Estimation Procedure.
# num_sims <- 100000
# dim <- 4
# all_results <- matrix(NA, nrow = num_sims, ncol = 4)
# norm_vcov <- 3.12 * diag(dim)
# norm_vcov / sqrt(diag(norm_vcov) %*% t(diag(norm_vcov)) )
# for (bb in 1:num_sims) {
#   if(bb %% 1000 == 0){cat(bb, " ")}
#   my_dat <- MASS::mvrnorm(n = 101, mu = rep(0, dim), Sigma = norm_vcov)
#   est <- sqrt(nrow(my_dat)) * est_pearson(my_dat)[1, 1]
#   ic_ests <- list(
#     "old" = est_influence_pearson(my_dat),
#     "new" = new_est_influence_pearson(my_dat),
#     "newer" = new_new_est_influence_pearson(my_dat)
#   )
#   var_ests <- lapply(ic_ests, FUN = function(dat) t(dat) %*% dat / nrow(dat))
#   var_ones <- sapply(var_ests, "[", 1, 1)
#   all_results[bb, ] <- c(est, var_ones)
# }
#
# apply(all_results, 2, mean, na.rm = TRUE)
# apply(all_results, 2, var, na.rm = TRUE)
#
# df_res <- as.data.frame(all_results)
# colnames(df_res) <- c("est", "var_est_old",
#                       "var_est_new", "var_est_newer")
#
# library(tidyverse)
# df_res <- df_res %>%
#   mutate("test_old" = 2 * (1 - pnorm(abs(est / sqrt(var_est_old)))),
#          "test_newer" = 2 * (1 - pnorm(abs(est / sqrt(var_est_newer)))))
# est_influence_pearson <- function(observ, trans = "none"){
# n <- nrow(observ)
# num_cov <- ncol(observ) - 1
# means <- colMeans(observ)
# y_mean <- means[1]
# cent_obs <- observ - matrix(rep(colMeans(observ), each = n), nrow = n)
# cent_obs_sqrd <- cent_obs^2
# sigmas <- colSums(cent_obs_sqrd)/(n - 1)
# covs <- as.numeric(crossprod(cent_obs[, 1], cent_obs[, -1]))/(n - 1)
# y_var <- sigmas[1]
# cent_cov <- cent_obs[, -1] * cent_obs[, 1] - matrix(rep(covs, each = n),
#                                                     nrow = n)
# cent_var <- cent_obs_sqrd - rep(sigmas, each = n)
# psi_1 <- matrix(rep(sigmas[-1], each = n), nrow = n) *
#   cent_var[, 1] + y_var * cent_var[, -1]
# ic <- (cent_cov - psi_1 * rep(covs/(y_var * sigmas[-1]), each = n)) *
#   rep((1/sqrt(y_var * sigmas[-1])), each = n)
# if(trans == "none"){
#   return(ic)
# }else if(trans == "tsqd"){
#   num_var <- ncol(observ)
#   rho <- stats::cor(observ[, 1], observ[, -1], method = "pearson")[1, ]
#   mat_mult <- diag(2 * rho / (1 - rho ** 2) ** 2)
#   return(ic %*% t(mat_mult))
# }
# }
#
# new_est_influence_pearson <- function(observ, trans = "none"){
#   n <- nrow(observ)
#   num_cov <- ncol(observ) - 1
#   means <- colMeans(observ)
#   y_mean <- means[1]
#   cent_obs <- observ - matrix(rep(colMeans(observ), each = n), nrow = n)
#   cent_obs_sqrd <- cent_obs ** 2
#   sigmas <- colSums(cent_obs_sqrd)/n
#   covs <- as.numeric(crossprod(cent_obs[, 1], cent_obs[, -1]))/n
#   var_x <- sigmas[-1]
#   y_var <- sigmas[1]
#   cent_cov <-
#     sweep(cent_obs[, -1], MARGIN = 1, cent_obs[, 1], `*`) -
#     matrix(rep(covs * (1 - 1 / n), each = n), nrow = n)
#   cent_var <- cent_obs_sqrd - rep(sigmas * (1 - 1 / n), each = n)
#   my_var <-
#     sweep(cent_cov, MARGIN = 2, (1 / sqrt(y_var * sigmas[-1])), `*`) -
#     sweep(cent_var[, -1], MARGIN = 2,
#           (1 / 2) * covs / sqrt((y_var * var_x^3)), `*`) -
#     sweep(matrix(rep(cent_var[, 1, drop = FALSE], num_cov), ncol = num_cov),
#           MARGIN = 2, (1 / 2) * covs / sqrt((y_var^3 * var_x)), `*`)
#   if(trans == "none"){
#     return(my_var)
#   }else if(trans == "tsqd"){
#     num_var <- ncol(observ)
#     rho <- stats::cor(observ[, 1], observ[, -1], method = "pearson")[1, ]
#     mat_mult <- diag(2 * rho / (1 - rho ** 2) ** 2)
#     return(ic %*% t(mat_mult))
#   }
# }
