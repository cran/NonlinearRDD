#' Standard error for nonlinear RDD
#'
#' Calculates the standard error for the semiparametric estimates in the nonlinear structural function in RDD
#'
#' @param df Data frame whose columns are outcome Y, treatment Treat, and running variable R, respectively.
#' @param gamma Value of the parameters in the structural function.
#' @param t Treatment value for calculating the marginal effect.
#' @param e Value of the error term for calculating the marginal effect.
#' @examples
#' \donttest{
#' n <- 500
#' set.seed(123)
#' data <- nrddgp(n)
#' nrdse(data,c(1,1,1),.5,.5)
#' }
#' @export

nrdse <- function(df,gamma,t,e){
  Y <- df[,1]
  Treat <- df[,2]
  R <- df[,3]
  n <- nrow(df)
  R0 <- R[R<0]
  R1 <- R[R>0]
  T0 <- Treat[R<0]
  T1 <- Treat[R>0]
  Y0 <- Y[R<0]
  Y1 <- Y[R>0]
  df.below <- df[df[,3]<0,]
  df.above <- df[df[,3]>0,]
  m1 <- lpcde::lpcde(T0,Y0)
  est1 <- m1$Estimate[,3]
  m0 <- lpcde::lpcde(T1,Y1)
  est0 <- m0$Estimate[,3]
  fR <- rddensity::rddensity(df[,3])$hat[[1]]
  se <- 10 *sum((est1 - est0)^2) / fR
  value <- c(1,t,e)
  M <- outer(value,value)
  se <- se * M
  return(se)
}
