#' This is an internal function. The kernel function
#'
#' @param x A numeric argument.
#' @return The value of the kernel function at x.

k <- function(x){
  return( (2 * abs(x)^3 - 3* x^2 + 1) * (abs(x) < 1) )
}
