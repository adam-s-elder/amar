#' A function used to calculate various L_p norms
#'
#' @export

l_p_norm <- function(x, p = "max", type = "lp"){
  if(type == "lp"){
    if (p == "max") {
      return(max(abs(x)))
    } else {
      l_p <- as.integer(p)
      return(sum(abs(x)**l_p)**(1 / l_p))
    }
  }else if (type == "ordl2"){
    l_p <- as.integer(p)
    some_x <- x[order(abs(x)), decreasing = TRUE]
    return(sum(some_x[1:p] **2))
  }
}
