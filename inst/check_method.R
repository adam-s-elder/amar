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
