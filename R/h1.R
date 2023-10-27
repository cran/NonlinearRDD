#' This is an internal function. The treatment choice function above the cutoff, used for the data generating function nrddgp.
#'
#' @param r Value of the running variable.
#' @param u Value of the conditional rank of the treatment given the running variable.
#' @return Value of the treatment choice function above the cutoff.

h1 <- function(r,u){
  return( r+2*u^3)
}
