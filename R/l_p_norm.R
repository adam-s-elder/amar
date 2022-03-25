#' A function used to calculate various L_p norms
#' @param x Observed data
#' @param p index of the norm
#' @param type Kind of norm used (currently only lp and
#' sum of squares norms are supported)
#' @return The norm of x, of type \code{type} of index p.  For example,
#' the euclidean norm x has \code{p = 2, type = "lp"}
#'
#' @examples
#' x <- c(3, 4)
#' l_p_norm(x, p = 2, type = "lp")
#' l_p_norm(x, p = 2, type = "ssq")
#' l_p_norm(x, p = "max", type = "lp")
#' l_p_norm(x, p = 1, type = "ssq")
#'
#' y <- c(3, 4, 5, 6)
#' l_p_norm(y, p = 4, type = "lp")
#' l_p_norm(y, p = 3, type = "ssq")
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
    return(sqrt(sum(some_x[1:p])))
  }
}
