#' This a function used to look at the IC.
#' 
#' @param IC an influence curve evaluated at each observation
#' @return Histogram of the IC characteristics
#'
#' @export


look_IC <- function(IC){
  colnames(IC) <- 1:ncol(IC)
  tidy_ic <- tidyr::pivot_longer(
    as.data.frame(IC), cols = colnames(IC),
    names_to = "Column")
  mean_df <- tidy_ic %>% group_by(Column) %>%
    summarise(mean = mean(value))
  
  tidy_ic %>%
    ggplot(aes(x = value)) + geom_histogram() + 
    geom_vline(data = mean_df, aes(xintercept = mean)) +
    facet_wrap(~as.numeric(Column))
}