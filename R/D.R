#' This is an internal function. It calculates the criterion function for semiparametric estimation.
#'
#' @param ga Parameter in the structural function.
#' @param Y Outcome variable.
#' @param Treat Treatment variable.
#' @param R Running variable.
#' @param b1 Bandwidth for the local linear estimation of the conditional distribution of Y given T,R.
#' @param b2 Bandwidth for the smoothing of the indicator function 1(Y <= y).
#' @param b3 Bandwidth for the first-step conditional quantile estimation of T given R.
#' @param Nonseparable Whether the structural function should be estimated as a separable or nonseparable function of the treatment. The default is TRUE.
#' @param e_space The space of epsilon on which the numerical integration is conducted.
#' @param u_grid Grid length of the quantile space for numerical integration. The default is 0.02.
#' @param t_tilde Normalization constant in the parametric specification of the structural function.
#' @return Returns the value of the criterion function for a given value of the parameter.

D <- function(ga,b1,b2,b3,Y,Treat,R,
              u_grid = 0.02,
              e_space,
              t_tilde,
              Nonseparable=TRUE)
{
  R0 <- R[R<0]
  R1 <- R[R>0]

  T0 <- Treat[R<0]
  T1 <- Treat[R>0]

  Y0 <- Y[R<0]
  Y1 <- Y[R>0]

  n <- length(R)
  n0 <- length(R0)
  n1 <- length(R1)

  e_seq <- e_space
  u_seq <- seq(0.05,0.95,by = u_grid)

  W0 <- k( R0 / b3 )
  W1 <- k( R1 / b3 )
  h0_hat <- quantreg::rq(T0 ~ 1, tau = u_seq, weights = W0)
  h1_hat <- quantreg::rq(T1 ~ 1, tau = u_seq, weights = W1)

  h0_est <- h0_hat$coefficients[1,]
  h1_est <- h1_hat$coefficients[1,]

  ll <- Vectorize(function(e,i){
    t0 <- h0_est[i]
    t1 <- h1_est[i]
    y0 <- g(ga,t0,e,t_tilde,Nonseparable)
    y1 <- g(ga,t1,e,t_tilde,Nonseparable)
    a <- Rfast::lmfit(x=cbind(rep(1,n0),T0-t0,R0),
               y=K_Y((y0-Y0)/b2),
               w=k((T0-t0)/b1)*k( R0 / b1 ) )$be[1]
    b <- Rfast::lmfit(x=cbind(rep(1,n1),T1-t1,R1),
               y=K_Y((y1-Y1)/b2),
               w=k((T1-t1)/b1)*k( R1 / b1 ) )$be[1]
    return(a-b)

  })

  D_eu <- outer(e_seq,1:length(u_seq),ll)
  D_eu_int <- apply(D_eu,1,cumsum)^2
  return(mean(t(D_eu_int) * stats::dnorm(e_seq)))

  rdrobust1 <- rdrobust::rdrobust(1,1)
  rddensity1 <- rddensity::rddensity(1)
  lpdensity1 <- lpdensity::lpdensity(1)
  lpcde1 <- lpcde::lpcde(1)
  rm(rdrobust1,rddensity1,lpdensity1,lpcde1)

}
