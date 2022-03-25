#' Function for calculating the influence function used for
#' the real data example.
#' @param obs_data the observed data.  The first column should be the outcome.
#' @param what the desired return value. Should be one of \code{"ic"}
#' (influence curve), \code{"est"} (estimate), or \code{"both"}.
#' @param control any other control parameters to be passed to the estimator.
#' @return If \code{what} is
#'
#' - \code{"est"}, then return the estimated parameter.
#'
#' - \code{"ic"}, then return the estimated IC of the parameter estimate.
#'
#' - \code{"both"}, then return both the parameter estimate and
#' corresponding estimated IC.
#'
#' @examples
#'
#' expit <- function(x) exp(x) / (1 + exp(x))
#' ws <- matrix(rnorm(3000), ncol = 3)
#' probs <- expit(ws  %*% c(-1, 0, 2))
#' y <- rbinom(n = nrow(probs), size = 1, prob = probs[, 1])
#' wts <-   abs(rnorm(length(y))) + 1
#' wts <- length(wts) * wts / sum(wts)
#' cats <- rep(1:10, 100)
#' obs_dat <- cbind(y, "cat" = cats, "wt" = wts, ws)
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

ic.data.examp <- function(obs_data, what = "both", control = NULL) {
  y_vc <- obs_data[, "y"]
  wts <- obs_data[, "wt"]
  cc_cat <- obs_data[, "cat"]
  cat_unique <- unique(cc_cat)
  ws <- obs_data[, -which(colnames(obs_data) %in% c("y", "wt", "cat"))]
  est <- rep(NA, ncol(ws))
  fin_IC <- matrix(NA, nrow = length(y_vc), ncol = ncol(ws))
  for (cov_idx in seq_len(ncol(ws))) {
    w_j <- ws[, cov_idx, drop = TRUE]
    glm_fit <- suppressWarnings(stats::glm(y_vc ~ w_j, weights = wts,
                   family = "binomial"))
    est[cov_idx] <- glm_fit$coefficients[2]
    beta_0 <- glm_fit$coefficients[1]
    beta_1 <- glm_fit$coefficients[2]
    m_b <- expit(beta_0 + beta_1 * w_j)
    m_b_prod <- -1 * m_b * (1 - m_b)
    des_mat <- cbind(1, w_j)
    mult_mat <- matrix(rep(m_b_prod * wts, times = 2),
                       ncol = 2, byrow = FALSE)
    m_n_inv <- solve(
      t(des_mat * mult_mat) %*% des_mat / 150
    )
    grad <- des_mat
    grad[, 1] <- grad[, 1] * (y_vc - m_b)
    grad[, 2] <- grad[, 2] * (y_vc - m_b)
    cross_prod <-  -1 * m_n_inv %*% t(grad)
    prod_term_b2 <- cross_prod[2, , drop = TRUE] #(a_j0 + b_j0)(y - m_beta_j0)
    part_one_beta_2 <- wts * prod_term_b2
    part_two_beta_2 <- rep(NA, length(part_one_beta_2))
    for (this_cat in cat_unique) {
      for (inf_stat in c(0, 1)){
        w_m <- which((cc_cat == this_cat) & (y_vc == inf_stat))
        part_two_beta_2[w_m] <- (1 - wts[w_m]) * mean(prod_term_b2[w_m])
      }
    }
    if (any(is.na(part_two_beta_2))) {warning(
      "Some of the case control categories were
      not covered in the calculation of the IC"
    )}
    fin_IC[, cov_idx] <- part_one_beta_2 + part_two_beta_2
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
