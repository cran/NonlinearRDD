#' This is an internal function. The anti-derivative of the kernel function.
#'
#' @param x A numeric argument.
#' @return The value of the anti-derivative function at x.

K_Y <- function(x)
{
  return( (x < 0 & x >-1) * (-x^4/2 -x^3 + x + 1/2) + (x >= 0 & x < 1) * (x^4/2 -x^3 + x + 1/2) + (x >= 1))
}
