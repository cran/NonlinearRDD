#' Nonlinear RDD estimator
#'
#' Semiparametric estimation of the possibly nonlinear and nonseparable structural function in regression discontinuity designs with a continuous treatment variable
#'
#' @param df Data frame whose columns are outcome Y, treatment Treat, and running variable R, respectively.
#' @param b1 Bandwidth for the local linear estimation of the conditional distribution of Y given T,R. Default is n^(-1/5), where n is the sample size.
#' @param b2 Bandwidth for the smoothing of the indicator function 1(Y <= y).
#' @param b3 Bandwidth for the first-step conditional quantile estimation of T given R.
#' @param t_tilde Normalization constant in the structural function.
#' @param e_space The space of epsilon on which the numerical integration is conducted.
#' @param Nonseparable Whether the structural function should be estimated as a separable or nonseparable function of the treatment. The default is TRUE.
#' @return The estimated parameters in the structural function.
#' @examples
#' n <- 20
#' set.seed(123)
#' Treat <- runif(n)
#' R <- runif(n)-0.5
#' epsilon <- Treat^2 + R
#' Y <- Treat + R + epsilon
#' df <- data.frame(Y=Y,Treat=Treat,R=R)
#' NRDD(df,e_space = seq(0,1,by=0.5),Nonseparable=FALSE)
#'
#' \donttest{
#' n <- 500
#' set.seed(123)
#' data <- nrddgp(n)
#' NRDD(data,e_space = seq(-10,10,by=1))
#' }
#' @export

NRDD <- function(df,b1=2*nrow(df)^{-1/5},
                 b2=2*nrow(df)^{-1/5},
                 b3=2*nrow(df)^{-1/5},
                 t_tilde=stats::median(df[,2]),
                 e_space,
                 Nonseparable=TRUE)
{
  Y <- df[,1]
  Treat <- df[,2]
  R <- df[,3]
  if (Nonseparable) {initial <- c(1,1,1)}
  else {initial <- c(1,1)}
  return(stats::nlm(D,initial,steptol=1e-5,gradtol=1e-5,
             b1=b1,
             b2=b2,
             b3=b3,
             Y=Y,
             Treat=Treat,
             R=R,
             e_space = e_space,
             t_tilde=t_tilde,
             Nonseparable=Nonseparable)$estimate)
}
