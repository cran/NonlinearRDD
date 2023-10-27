#' Data generating process for nonlinear RDD
#'
#' Generates a sample of observations based on RDD with a continuous treatment variable and a possibly nonlinear and nonseparable structural function.
#'
#' @param n Sample size.
#' @param rho_R Correlation between R and epsilon.
#' @param rho_U Correlation between U and epsilon.
#' @param gamma_star true parameter value for the parameter in the structural function.
#' @param Nonseparable Whether the structural function should be estimated as a separable or nonseparable function of the treatment. The default is TRUE.
#' @return Data frame whose columns are outcome Y, treatment Treat, and running variable R, respectively.
#' @examples
#' n <- 500
#' set.seed(123)
#' data <- nrddgp(n)
#' @export

nrddgp <- function(n,rho_R=0.3,rho_U=0.3,gamma_star=c(1,1,1),Nonseparable=TRUE){
  dist_ReU <- copula::mvdc(copula::normalCopula(param = c(rho_R,0,rho_U),dim=3,dispstr = "un"),margins = c("norm","beta","unif"),
                   paramMargins = list(list(mean=0,sd=1),list(shape1=2,shape2=2),list(min=0,max=1)))
  ReU <- copula::rMvdc(n,dist_ReU)
  R <- ReU[,1]
  R0 <- R[R<0]
  R1 <- R[R>0]
  e0 <- ReU[R<0,2]
  e1 <- ReU[R>0,2]
  U0 <- ReU[R<0,3]
  U1 <- ReU[R>0,3]

  T0 <- h0(R0,U0)
  T1 <- h1(R1,U1)
  t_tilde <- 0.5

  Y0 <- g(gamma_star,T0,e0,t_tilde) + R0
  Y1 <- g(gamma_star,T1,e1,t_tilde) + R1

  Y <- c(Y0,Y1)
  Treat <- c(T0,T1)
  df <- data.frame(Y=Y,Treat=Treat,R=R)
  return(df)
}
