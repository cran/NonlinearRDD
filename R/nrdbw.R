#' Bandwidth selector for nonlinear RDD
#'
#' Provides recommendations for the three bandwidths in the semiparametric estimation of the structural function in RDD with a continuous treatment variable.
#'
#' @param df Data frame whose columns are outcome Y, treatment Treat, and running variable R, respectively.
#' @return A vector of length 3 containing the recommended bandwidth for each steps in the NRDD function.
#' @examples
#' n <- 500
#' set.seed(123)
#' data <- nrddgp(n)
#' nrdbw(data)
#' @export

nrdbw <- function(df){
  Y <- df[,1]
  Treat <- df[,2]
  R <- df[,3]
  n <- nrow(df)

  b1 <- rdrobust::rdbwselect(as.numeric(Y <= stats::median(Y)),R)$bws[1] * n^{1/5 - 1/4.5}
  b2 <- b1
  b3mean <- rdrobust::rdbwselect(Treat,R)$bws[1]
  b3 <- b3mean * (0.5^2 * 2 *pi)^{1/5}
  return(c(b1,b2,b3))
}
