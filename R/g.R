#' This is an internal function. The structural function.
#'
#' @param ga Value of the parameter.
#' @param t Treatment value.
#' @param e Unobserved heterogeneity value.
#' @param t_tilde Normalization constant.
#' @param Nonseparable Whether the structural function should be estimated as a separable or nonseparable function of the treatment. The default is TRUE.
#' @return The value of the structural function.

g <- function(ga,t,e,t_tilde,Nonseparable=TRUE){
  if (Nonseparable) {y <- ga[1] * (t-0.5) + ga[2] * (t^2 - t_tilde^2) + ga[3] * (t - t_tilde) * e + e}
  else {y <- ga[1] * (t-0.5) + ga[2] * (t^2 - t_tilde^2) + e}
  return(y)
}
