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
#' calculated.
#' @param observ the observed data
#' @param trans If there is a transformation of the parameter of interest
#'
#' @export

est_influence_pearson <- function(observ, trans = "none"){
  n <- nrow(observ)
  num_cov <- ncol(observ) - 1
  means <- colMeans(observ)
  y_mean <- means[1]
  cent_obs <- observ - matrix(rep(colMeans(observ), each = n), nrow = n)
  cent_obs_sqrd <- cent_obs^2
  sigmas <- colSums(cent_obs_sqrd)/(n - 1)
  covs <- as.numeric(crossprod(cent_obs[, 1], cent_obs[, -1]))/(n - 1)
  y_var <- sigmas[1]
  cent_cov <- cent_obs[, -1] * cent_obs[, 1] - rep(covs, each = n)
  cent_var <- cent_obs_sqrd - rep(sigmas, each = n)
  psi_1 <- matrix(rep(sigmas[-1], each = n), nrow = n) *
    cent_var[, 1] + y_var * cent_var[, -1]
  ic <- (cent_cov - psi_1 * rep(covs/(y_var * sigmas[-1]), each = n)) *
    rep((1/sqrt(y_var * sigmas[-1])), each = n)
  if(trans == "none"){
    return(ic)
  }else if(trans == "tsqd"){
    num_var <- ncol(observ)
    rho <- stats::cor(observ[, 1], observ[, -1], method = "pearson")[1, ]
    mat_mult <- diag(2 * rho / (1 - rho ** 2) ** 2)
    return(ic %*% t(mat_mult))
  }
}

est_pearson <- function(observ, trans = "none"){
    num_var <- ncol(observ)
    rho <- stats::cor(observ[, 1], observ[, -1], method = "pearson")
    if(trans == "tsqd"){
      return(rho ** 2 / (1 - rho ** 2))
    }else{
      return(rho)
    }
}

pearson <- function(est_or_IC){
  if(est_or_IC == "est"){return(est_pearson)}
  if(est_or_IC == "IC"){return(est_influence_pearson)}
  else{stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}

#### Placing all of the different influence function and parameter estimation
#### inside of a single function

#' Wrapper function that provides all of the possible influence functions
#'
#' @param param name of the parameter of interest
#' @param est_or_IC specify if an estimator of the parameter or IC is needed
#'
#' @export


get_infl <- function(param, est_or_IC){
  if(param == "spearman"){return(spearman(est_or_IC))}
  if(param == "pearson"){return(pearson(est_or_IC))}
  if(param == "mean_f"){return(mean_f(est_or_IC))}
  if(param == "median"){return(median(est_or_IC))}
  if(param == "DE2"){return(de_2(est_or_IC))}
  if(param == "DE2nl"){return(de_2_nl(est_or_IC))}
  if(param == "DE2KM"){return(de_2_know_mis(est_or_IC))}
  if(param == "DE3"){return(de_3(est_or_IC))}
  if(param == "RDE"){return(rde(est_or_IC))}
  else{
    stop("The specified parameter is not one of the included parameters,
         please choose from 'spearman', 'pearson', 'mean_f', 'DE2', or 'DE3'")
  }
}
