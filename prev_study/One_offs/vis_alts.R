## This script was written in order to look at what detectable alternatives looked like
## visually (the distributions look very similar).  I was also interested in looking at
## How similar the histograms of the distance from the null and alternative look between
## the two generated datasets (null and alternative).  

null_dstr <- MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = matrix(c(1, 0.3, 1, 0.3), nrow = 2))
alt_dstr <- MASS::mvrnorm(n = 100, mu = rep(0.05, 2), Sigma = matrix(c(1, 0.3, 1, 0.3), nrow = 2))
new_obs <- function(x) return(c((x[1] ** 2 + x[2] ** 2)**(1/2), 
                                ((x[1] - 0.05) ** 2 + (x[2] - 0.05) ** 2)**(1/2)))
shft_null <- t(apply(null_dstr, 1, new_obs))
shft_alt <- t(apply(alt_dstr, 1, new_obs))
all_obs <- rbind(shft_null, shft_alt)
plot(all_obs, col = rep(c("blue", "red"), each = 100))
all_data <- as.data.frame(cbind("Y" = rep(c(0, 1), each = 100), all_obs))
summary(lm(Y~., data = all_data))
svm_obj <- e1071::svm(Y~., data = all_data)
lines(c(0, 3), c(-0.025, 3- 0.025))
plot(rbind(null_dstr, alt_dstr), col = rep(c("blue", "red"), each = 100))


null_dstr <- MASS::mvrnorm(n = 100000, mu = rep(0, 2), Sigma = matrix(c(1, 0.3, 1, 0.3), nrow = 2))
alt_dstr <- MASS::mvrnorm(n = 100000, mu = rep(0.05, 2), Sigma = matrix(c(1, 0.3, 1, 0.3), nrow = 2))
shft_null <- t(apply(null_dstr, 1, new_obs))
shft_alt <- t(apply(alt_dstr, 1, new_obs))
apply(shft_null, 2, mean)
apply(shft_alt, 2, mean)
mean(shft_null[, 1] > shft_null[, 2])
mean(shft_alt[, 1] > shft_alt[, 2])

mean()