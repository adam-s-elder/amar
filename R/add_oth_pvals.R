#' Add Other p-values
#'
#' Add pvalues for the Liu and Xie and Bonferroni based
#' tests using the estimated parameter estimates and
#' corresponding standard errors.
#'
#' @param test_result The test result from mv_pn_test
#'
#' @return The same test result object with additional p-values for
#' the Liu and Xie (2021) test (liu_xie_pvalue) and the Bonferroni based
#' test (bonf_pvalue)
#' @export
#'
#' @examples
#' ## NOTE: More monte-carlo samples should be taken are taken here.  This is
#' ## only done to lower computation time.
#' set.seed(10)
#' test <- mv_pn_test(data.frame(y = rnorm(100), x = rnorm(100)),
#'                    ic.pearson, test.control(n_peld_mc_samples = 20,
#'                                             ts_ld_bs_samp = 20))
#' test_with_extra <- add_oth_pvals(test)
#' test_with_extra[c("pvalue", "liu_xie_pvalue", "bonf_pvalue")]

add_oth_pvals <- function(test_result) {
  new_results <- test_result
  zscores <- abs(test_result$param_ests / test_result$param_ses)
  # Liu and Xie p-value
  cauchy_pvalue <-  1 - stats::pcauchy(mean(tan((
    2 * stats::pnorm(zscores) - 3/2) * pi
    )))
  new_results$liu_xie_pvalue <- cauchy_pvalue
  # Bonferroni based p-value
  bonf <- min(2 * stats::pnorm(-zscores)) * length(zscores)
  new_results$bonf_pvalue <- bonf
  return(new_results)
}
