#' Function for calculating the influence function used for
#' the real data example.
#' @param obs_data the observed data.  The first column should be the outcome.
#' @param what the desired return value. Should be one of `"ic"`
#' (infludence curve), `"est"` (estimate), or `"both"`.
#' @param control any other control parameters to be passed to the estimator.
#'
#' @examples
#'
#' expit <- function(x) exp(x) / (1 + exp(x))
#' ws <- matrix(rnorm(30000), ncol = 3)
#' probs <- expit(ws  %*% c(-1, 0, 2))
#' y <- rbinom(n = nrow(probs), size = 1, prob = probs[, 1])
#' wts <-   abs(rnorm(length(y))) + 1
#' wts <- length(wts) * wts / sum(wts)
#' obs_dat <- cbind(y, "wt" = wts, ws)
#' est_ic <- ic.data.examp(obs_dat, what = "both")
#' my_est <- est_ic$est
#' my_ic <- est_ic$ic / nrow(ws)
#' var_mat <- t(my_ic) %*% my_ic
#' sqrt(diag(var_mat))
#' for(cov_idx in 1:ncol(ws)){
#'  print(summary(stats::glm(y ~ ws[, cov_idx], weights = obs_dat[, "wt"],
#'                     family = binomial))$coefficients[2, 1:2])
#' }
#'
#' @export

ic.data.examp <- function (obs_data, what = "both", control = NULL) {
  y_vc <- obs_data[, "y"]
  wts <- obs_data[, "wt"]
  ws <- obs_data[, - which(colnames(obs_data) %in% c("y", "wt"))]
  est <- rep(NA, ncol(ws))
  fin_IC <- matrix(NA, nrow = length(y_vc), ncol = ncol(ws))

  for (cov_idx in seq_len(ncol(ws))) {
    w_j <- ws[, cov_idx]
    glm_fit <- suppressWarnings(stats::glm(y_vc ~ w_j, weights = wts,
                   family = "binomial"))
    est[cov_idx] <- glm_fit$coefficients[2]
    beta_0 <- glm_fit$coefficients[1]
    beta_1 <- glm_fit$coefficients[2]
    jst_expit <- expit(beta_0 + beta_1 * w_j)
    expit_fracs <-
      jst_expit / (
        1 + exp(beta_0 + beta_1 * w_j)
      )

    exp_w2 <- mean(expit_fracs * (w_j ** 2) * wts) #I_{beta,11}
    exp_w <- mean(expit_fracs * w_j * wts) #I_{beta,21} and I_{beta,21}
    exp_1 <- mean(expit_fracs * wts) # I_{beta, 22}

    c_b <- (exp_w2 * exp_1 - (exp_w ** 2)) ** (-1)
    ic_vals <- c_b * wts * (y_vc - jst_expit) * (w_j * exp_1 - exp_w)
    fin_IC[, cov_idx] <- ic_vals
  }
  ret <- list()
  if (what %in% c("both", "est")) {
    ret$est <- est
  }
  if (what %in% c("both", "ic")) {
    ret$ic <- fin_IC
  }
  return(ret)
}

expit <- function(x) exp(x) / (1 + exp(x))
