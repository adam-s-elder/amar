#' Generate data using one of the four specified models
#' from MCKEAGUE and QIAN paper
#' @param ss sample size (number of generated observations)
#' @param dim dimension of the generated data
#' @param rho between x correlation
#' @param model number indicating which model should be used
#' @param b local alternative to be used
#' @return A dataframe with outcome in the first column and covariates
#' in the others.  Data are generated according to the models described
#' originally in McKeague and Qian (2015).
#' @examples
#' # Sample data
#' null_dat <- make_data(ss = 100, dim = 10, rho = 0, model = 1)
#' model_two_dat <- make_data(ss = 100, dim = 10, rho = 0, model = 2)
#' model_three_dat <- make_data(ss = 100, dim = 10, rho = 0, model = 2)
#'
#' @export

make_data <- function(ss, dim, rho, model = 1, b = NULL){
  x_cov <- matrix(rho, nrow = dim, ncol = dim)
  diag(x_cov) <- 1
  X <- MASS::mvrnorm(n = ss, mu = rep(0, dim), Sigma = x_cov)
  epsilon <- stats::rnorm(ss)
  if (model == 1) {
    data <- cbind(epsilon, X)
  }
  if (model == 2) {
    data <- cbind(X[, 1]/4 + epsilon, X)
  }
  if (model == 3) {
    beta <- c(rep(c(0.15, -0.1), each = 5), rep(0, dim - 10))
    data <- cbind(X %*% beta + epsilon, X)
  }
  if (model == 4) {
    # A potential other alternative to consider
    # though this is not considered anywhere before.
    theta_n <- c(b, rep(0, dim - 1)) / sqrt(ss)
    data <- cbind(X %*% theta_n + epsilon, X)
  }
  return(data)
}
