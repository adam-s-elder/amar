#' A helper function to generate a multiplier bootstrap sample
#'
#' @param n Number of desired observations from your bootstrapped sample.
#' @param obs_ic The empirical estimate of the influence curve to be
#' used in the multiplier bootstrap.
#' @param param_est The estimated parameters which will be the mean of the
#' multiplier bootstrap sample.
#' @param epsilon_mat The matrix of Normal observations with independent
#' observations from a normal with an identity Covariance matrix.
#' @param center Boolean.  If true, the bootstrapped data will be centered at
#' zero.  Otherwise, it will be centered at param_est.
#' @param rate Normalizing constant. Should either be \code{"n"} or
#' \code{"rootn"}.
#' @return A sample of size \code{n} generated using a multiplier bootstrap
#' with a variance given by t(\code{obs_ic})%*%\code{obs_ic}.
#'
#' @export

gen_boot_sample <- function(epsilon_mat , obs_ic, center = TRUE,
                            param_est = 0, rate = "n"){
  num_obs <- nrow(obs_ic)
  if (rate == "rootn") {
    cent_boot <-  epsilon_mat %*% obs_ic / sqrt(num_obs)
  }else if (rate == "n") {
    cent_boot <-  epsilon_mat %*% obs_ic / num_obs
  }
  if (center == TRUE) {
    return(cent_boot)
  }else{
    non_cent_boot <- sweep(cent_boot, 1, param_est, "+")
    return(non_cent_boot)
  }
}
