#' Control function for the adaptive norm test
#'
#' @param big_train Data is split into splits of roughly equal sizes. The
#' number of splits is equal to num_folds. If big_train is TRUE then all
#' but one of these splits will be training data, if big_train is
#' FALSE all but one will be testing data.
#' @param f_cv_summary How test statistics from different folds are combined to create an overall
#' test statistic.  Usually the mean is used.
#' @param more_info Boolean indicating if the test schould return more
#' information that just the p-value.  When true, the chosen norm index,
#' the bonferroni based p-value, the test statistic,
#' and the estimated distribution of the test statistic will be returned.
#' @param n_bs_smp Number of samples to be used in estimating the limiting distribution of the
#' test statistic under the null.
#' @param nrm_type The type of norm to be used for the test.
#' Generally the l_p norm
#' @param nrmlize Boolean for if the estimator should be normalized to have
#' a limiting distribution with identity covaraince matrix.
#' @param num_folds The number of folds to be used in the cross-validation procedure.  If set to 1,
#' no cross validation will be used.
#' @param perf_meas the prefered measure used to generate the test statistic.
#' @param pos_lp_norms The index of the norms to be considered.  For example if we use the l_p norm,
#' norms_indx specifies the different p's to try.
#' @param show_hist set to True if you would like to see a histogram of the test statistic's
#' limiting distribution bootstrap draws compared to the estimated parameter.
#' @param test_stat_func A function that will provide the test statistic
#' for the given fold (using the testing data), and uses the best
#' norm (decided on using the training data).
#' @param ts_ld_bs_samp The number of test statistic limiting distribution bootstrap samples to be
#' drawn.
#' @param ld_est_meth Method for estimating the limiting distribution of the test statistic
#' parametric bootstrap or permutation (permutation works better for small sample sizes)
#' @param ... Other arguments needed in other places.
#'
#' @export

test.control <- function(
  big_train = TRUE,
  f_cv_summary = mean,
  more_info = NULL,
  n_bs_smp = 300,
  nrm_type = "lp",
  nrmlize = FALSE,
  num_folds = 1,
  perf_meas = "est_pow",
  pos_lp_norms = c(1, 2, 3, "max"),
  show_hist = FALSE,
  test_stat_func = "mag",
  ld_est_meth = "par_boot",
  ts_ld_bs_samp = 250, ...
){
  formal_args <- formals(sys.function())
  dot_args <- list(...)
  p <- .get.args(formal_args, dot_args)
  p$train_mlt <- (-1) ** (2 + as.integer(p$big_train))
  return(p)
}


## Pulled from epiModel
.get.args <- function (formal.args, dot.args)
{
  p <- list()
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg, pos = parent.frame()))
  }
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }
  return(p)
}
