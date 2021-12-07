#' Functions for creating data
#' @param n number of observations
#' @param d number of dimentions
#' @param cov desired covaraince
#' @export

gen_exp_cop <- function(n, d, cov){
  rate_params <- c(0.5, 1, 2)
  gendat <- cop_data <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = cov)
  unifdat <- stats::pnorm(gendat)
  for(i in 1:d){
    cop_data[, i] <- stats::qexp(unifdat[, i], rate = rate_params[1 + (i %% 3)])
  }
  return(cop_data)
}

#' Generate observations from a normal distribution
#' @param n number of observations
#' @param d number of covariates
#' @param cov variance-covariance of the multivariate normal distribution.
#' @export
#'

gen_norm <- function(n, d, cov) {MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = cov)}


#' Generate data from a poission random variable while still obtaining the desired
#' xy correlation for each xy match and obtain the desired xx correlation among all
#' x variables. The y and non-correlated x will be marginally
#' independent.
#'
#' @param base Base rate that is desired
#' @param num_cor Number of correlated X's
#' @param xycor the correlation between each correlated x and y
#' @param xxcor the correlation between the x's
#' @param num_uncor Number of uncorrelated X's
#' @export

gen_one_samp_mi <- function(base, xycor, xxcor, num_cor, num_uncor){
  zv <- stats::rpois(1, base)
  cr_x <- stats::rpois(num_cor, base)
  uncr_x <- stats::rpois(num_uncor, base)
  y <- base ; z <- xxcor
  a <- num_cor * xycor ** 2 - 1
  b <- (y * (num_cor + 1) + num_cor * z) * xycor ** 2
  c <- (y**2 + y * z) * xycor ** 2
  xy_lam <- (-b - sqrt(b**2 - 4 * a * c))/(2 * a)
  # xy_lam <- -xycor**2 * ((num_cor + 1) * y + num_cor * z) - sqrt((xycor**2 * (2 * y + z)) ** 2 - 4 * xycor **2 * num_cor * (xycor ** 2 - 1) * (y**2 + y * z))/(2 * num_cor * (xycor **2 - 1))
  # xy_lam <- base * xycor/(1 - xycor)
  xx_lam <- base * xxcor/(1 - xxcor)
  xy_add <- stats::rpois(num_cor, xy_lam)
  xx_add <- stats::rpois(1, xx_lam)
  zv <- zv + sum(xy_add)
  cr_x <- cr_x + xy_add + xx_add
  uncr_x <- uncr_x + xx_add
  return(c(zv, cr_x, uncr_x))
}

#' Generate data from a poission random variable while still obtaining the desired
#' xy correlation for each xy match and obtain the desired xx correlation among all
#' x variables. The y and x will always be marginally dependent due to the xx correlation.
#'
#' @param base Base rate that is desired
#' @param num_cor Number of correlated X's
#' @param xycor the correlation between each correlated x and y
#' @param xxcor the correlation between the x's
#' @param num_uncor Number of uncorrelated X's
#' @export

gen_one_samp_md <- function(base, xycor, xxcor, num_cor, num_uncor){
  uncr_x <- stats::rpois(num_uncor, base)
  xx_lam <- base * xxcor/(1 - xxcor)
  c1 <- num_cor * xx_lam + base
  c2 <- num_cor * (base + xx_lam) ** 2
  c3 <- (base + xx_lam) * num_cor * c1
  xy_lam <- (c2 * xycor ** 2) / (c1 ** 2 - c3 * xycor ** 2)
  # xy_lam_2 <- ((base + xx_lam) + (num_cor - 1) * (xx_lam))/(num_cor * xycor ** 2)
  cr_x <- stats::rpois(num_cor, base)
  xx_add <- stats::rpois(1, xx_lam)
  cr_x <- cr_x + xx_add
  uncr_x <- uncr_x + xx_add
  if(xy_lam != 0){
    zv <- stats::rpois(1, xy_lam * sum(cr_x))
  }else{
    zv <- stats::rpois(1, 10)
  }
  return(c(zv, cr_x, uncr_x))
}

#' Generate data from a poission random variable while still obtaining the desired
#' xy correlation for each xy match and obtain the desired xx correlation among all
#' x variables
#'
#' @param ss sample size.
#' @param base Base rate that is desired
#' @param xycor the correlation between each correlated x and y
#' @param xxcor the correlation between the x's
#' @param num_cor Number of correlated X's
#' @param num_uncor Number of uncorrelated X's
#' @param marg_indep Are the X-values independent from one another
#' @export

gen_pois_samp <- function(ss, base, xycor, xxcor, num_cor, num_uncor, marg_indep = TRUE){
  all_obs <- matrix(NA, nrow = ss, ncol = num_cor + num_uncor + 1)
  if (marg_indep){
    for(samp_idx in 1:ss){
      all_obs[samp_idx, ] <- gen_one_samp_mi(base, xycor, xxcor, num_cor, num_uncor)
    }
    return(all_obs)
  }else{
    for(samp_idx in 1:ss){
      all_obs[samp_idx, ] <- gen_one_samp_md(base, xycor, xxcor, num_cor, num_uncor)
    }
    return(all_obs)
  }
}


#' Generate data using one of the four specified models
#' from MCKEAGUE and QIAN paper
#' @param ss sample size (number of generated observations)
#' @param dim dimension of the generated data
#' @param rho between x correlation
#' @param model number indicating which model should be used
#' @param b local alternative to be used
#' @export
#'

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
    theta_n <- c(b, rep(0, dim - 1))/sqrt(ss)
    data <- cbind(X %*% theta_n + epsilon, X)
  }
  return(data)
}
