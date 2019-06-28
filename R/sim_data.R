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


gen_one_samp_mi <- function(base, xycor, xxcor, num_cor, num_uncor){
  zv <- rpois(1, base)
  cr_x <- rpois(num_cor, base)
  uncr_x <- rpois(num_uncor, base)
  y <- base ; z <- xxcor
  a <- num_cor * xycor ** 2 - 1
  b <- (y * (num_cor + 1) + num_cor * z) * xycor ** 2
  c <- (y**2 + y * z) * xycor ** 2
  xy_lam <- (-b - sqrt(b**2 - 4 * a * c))/(2 * a)
  # xy_lam <- -xycor**2 * ((num_cor + 1) * y + num_cor * z) - sqrt((xycor**2 * (2 * y + z)) ** 2 - 4 * xycor **2 * num_cor * (xycor ** 2 - 1) * (y**2 + y * z))/(2 * num_cor * (xycor **2 - 1))
  # xy_lam <- base * xycor/(1 - xycor)
  xx_lam <- base * xxcor/(1 - xxcor)
  xy_add <- rpois(num_cor, xy_lam)
  xx_add <- rpois(1, xx_lam)
  zv <- zv + sum(xy_add)
  cr_x <- cr_x + xy_add + xx_add
  uncr_x <- uncr_x + xx_add
  return(c(zv, cr_x, uncr_x))
}

gen_one_samp_md <- function(base, xycor, xxcor, num_cor, num_uncor){
  uncr_x <- rpois(num_uncor, base)
  xx_lam <- base * xxcor/(1 - xxcor)
  c1 <- num_cor * xx_lam + base
  c2 <- num_cor * (base + xx_lam) ** 2
  c3 <- (base + xx_lam) * num_cor * c1
  xy_lam <- (c2 * xycor ** 2) / (c1 ** 2 - c3 * xycor ** 2)
  # xy_lam_2 <- ((base + xx_lam) + (num_cor - 1) * (xx_lam))/(num_cor * xycor ** 2)
  cr_x <- rpois(num_cor, base)
  xx_add <- rpois(1, xx_lam)
  cr_x <- cr_x + xx_add
  uncr_x <- uncr_x + xx_add
  if(xy_lam != 0){
    zv <- rpois(1, xy_lam * sum(cr_x))
  }else{
    zv <- rpois(1, 10)
  }
  return(c(zv, cr_x, uncr_x))
}

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
