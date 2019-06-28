##################################################
#### Generating observations
####
#### This script is used to create observations
#### from various multivariate distributions
#### to test our method.
####
#### Author : Adam Elder
##################################################

gen_exp_cop <- function(n, d, cov){
  rate_params <- c(0.5, 1, 2)
  gendat <- cop_data <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = cov)
  unifdat <- pnorm(gendat)
  for(i in 1:d){
    cop_data[, i] <- qexp(unifdat[, i], rate = rate_params[1 + (i %% 3)])
  }
  return(cop_data)
}

gen_norm <- function(n, d, cov) {MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = cov)}


