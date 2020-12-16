#' @title An implementation of approximately compute the Jacobian Matrix.
#' @description This is an implementation of approximately compute the Jacobian
#'     Matrix of a multivariate function.
#' @param parvec a vector (x1, x2, ..., xk)(vector)
#' @param infcn the input function (g1, g2, ..., gm) (vector)
#' @param eps a small difference in parvec, with default value 1e-06 (numeric)
#' @return a jacobian matrix (gradient matrix) of function (g1, g2, ..., gm)
#' @examples
#' \dontrun{
#' f=function(x) c(x[1]^3-x[2]-1,x[1]^3-x[2]^2-1)
#' f(c(1,2))
#' Gradmat(c(1,2),f)
#' }
#' @export
Gradmat<-function(parvec, infcn, eps = 1e-06)
{
  dd = length(parvec) 
  aa = length(infcn(parvec)) 
  epsmat = (diag(dd) * eps)/2 
  gmat = array(0, dim = c(aa, dd))
  for(i in 1:dd)
    gmat[, i] = (infcn(parvec + epsmat[, i]) - 
                   infcn(parvec - epsmat[, i]))/eps 
  if(aa > 1) gmat else c(gmat) 
}

#' @title An implementation of Newton-Raphson Method for root finding for multivariate function.
#' @description An implementation of Newton-Raphson Method for root finding for multivariate function.
#' @param inipar an initial value for vector (x1_0, x2_0, ..., xk_0)(vector)
#' @param infcn the input function (g1, g2, ..., gm) (vector)
#' @param nmax the maximum iteration times, with default 25 (int)
#' @param stoptol the tolerance condition for stopping the iteration (numeric)
#' @param eps a small difference in parvec, with default value 1e-06 (numeric)
#' @param gradfunc the gradient matrix passed from outside, with default NULL (logical)
#' @return a list that contains: the iteration times, the root founded by the method, the initial condition, the function value at root.
#' @examples
#' \dontrun{
#' f=function(x) c(x[1]^3-x[2]-1,x[1]^3-x[2]^2-1)
#' NRroot(c(1,2),f)
#' }
#' @export
NRroot<-function(inipar, infcn, nmax = 50, stoptol = 1e-05,
         eps = 1e-06, gradfunc = NULL)
{
  if(is.null(gradfunc)) 
    gradfunc = function(x) Gradmat(x, infcn, eps)
  ctr = 0 
  newpar = inipar 
  oldpar = inipar - 1 
  while(ctr < nmax & sqrt(sum((newpar - oldpar)^2)) > stoptol) {
    oldpar = newpar
    newpar = oldpar - solve(gradfunc(oldpar), infcn(oldpar)) 
    ctr = ctr + 1
  }
  list(nstep = ctr, initial = inipar, final = newpar,
       funcval = infcn(newpar))
}

