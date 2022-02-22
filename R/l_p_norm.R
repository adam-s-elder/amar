#' A function used to calculate various L_p norms
#' @param x Observed data
#' @param p index of the norm
#' @param type Kind of norm used
#'
#' @export

l_p_norm <- function(x, p = "max", type = "lp"){
  if (!(type %in% c("lp", "ssq"))) {
    stop(paste0(
      "Currently the l_p_norm function only supports two types of norms, ",
      "including lp and ssq.  The norm type provided was ", type, "."
    ))
  }
  if (type == "lp"){
    if (p == "max") {
      return(max(abs(x)))
    } else {
      l_p <- as.integer(p)
      return(sum(abs(x)**l_p)**(1 / l_p))
    }
  }else if (type == "ssq") {
    l_p <- as.integer(p)
    x2 <- x ** 2
    some_x <- sort(x2, decreasing = TRUE)
    return(sum(some_x[1:p]))
  }
}
