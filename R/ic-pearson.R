#' Pearson Correlation IC and estimate
#'
#' This function takes a set of observations, and returns an estimate
#' and its corresponding estimated IC matrix for the estimates of the
#' pearson correlation. Estimates of the covariance are generated using
#' the empirical influence function.  The first column of your data should
#' correspond to the variable of interest (the variable for which pearson
#' correlation is calculated).
#' @param observ the observed data.  The first column should be the outcome.
#' @param what the desired return value. Should be one of \code{"ic"}
#' (influence curve), \code{"est"} (estimate), or \code{"both"}.
#' @param control any other control parameters to be passed to the estimator.
#' @return If \code{what} is:
#'
#' - \code{"est"}, then return the estimated person correlation.
#'
#' - \code{"ic"}, then return the estimated IC of the person correlation estimate.
#'
#' - \code{"both"}, then return both the estimated pearson correlation and the
#' estimated IC of the person correlation estimate.
#'
#' @examples
#' dat <- matrix(rnorm(80), nrow = 20)
#' ic.pearson(dat, what = "both")
#' ## Note that the estimate is the same as what is found using \code{cor}
#' cor(dat)[1, ]
#'
#' @export

ic.pearson <- function(observ, what = "both", control = NULL) {
  if (!(what %in% c("ic", "est", "both"))) {
    stop("what must be one of ic (influence curve), est (estimate), or both")
  }
  if (is.data.frame(observ)) {
    observ <- as.matrix(observ)
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
                                 obs_cent[, -1, drop = FALSE])) / n
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
  if (what %in% c("both", "est")) {
    ret$est <- stats::cor(observ[, 1], observ[, -1, drop = FALSE],
                          method = "pearson")
  }
  return(ret)
}
