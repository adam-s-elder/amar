#' This a function used to look at the IC.
#'
#' @param IC an influence curve evaluated at each observation
#' @return Histogram of the IC at each obeservation for each covariate.
#' Vertical lines indicate the mean of the IC.
#'
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(20)
#' ic <- amp::ic.pearson(matrix(rnorm(120), ncol = 6))$ic
#' look_IC(ic)
#'
#' @export

look_IC <- function(IC) {
  colnames(IC) <- 1:ncol(IC)
  tidy_ic <- tidyr::pivot_longer(
    as.data.frame(IC), cols = colnames(IC),
    names_to = "Column")
  mean_df <- dplyr::summarise(
    dplyr::group_by(tidy_ic, .data$Column), mean = mean(.data$value)
    )

  ggplot2::ggplot(tidy_ic, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram() +
    ggplot2::geom_vline(data = mean_df, ggplot2::aes(xintercept = .data$mean)) +
    ggplot2::facet_wrap(~as.numeric(.data$Column))
}
