## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp20065)

## ---- eval=FALSE--------------------------------------------------------------
#  Gradmat<-function(parvec, infcn, eps = 1e-06)
#  {
#    dd = length(parvec)
#    aa = length(infcn(parvec))
#    epsmat = (diag(dd) * eps)/2
#    gmat = array(0, dim = c(aa, dd))
#    for(i in 1:dd)
#      gmat[, i] = (infcn(parvec + epsmat[, i]) -
#                     infcn(parvec - epsmat[, i]))/eps
#    if(aa > 1) gmat else c(gmat)
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  NRroot<-function(inipar, infcn, nmax = 50, stoptol = 1e-05,
#           eps = 1e-06, gradfunc = NULL)
#  {
#    if(is.null(gradfunc))
#      gradfunc = function(x) Gradmat(x, infcn, eps)
#    ctr = 0
#    newpar = inipar
#    oldpar = inipar - 1
#    while(ctr < nmax & sqrt(sum((newpar - oldpar)^2)) > stoptol) {
#      oldpar = newpar
#      newpar = oldpar - solve(gradfunc(oldpar), infcn(oldpar))
#      ctr = ctr + 1
#    }
#    list(nstep = ctr, initial = inipar, final = newpar,
#         funcval = infcn(newpar))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix Gradmat_cpp(NumericVector parvec, Rcpp::Function infcn, double eps=1e-06){
#    int dd=parvec.size();
#    NumericVector temp=infcn(parvec);
#    int aa=(temp).size();
#    NumericMatrix epsmat;
#    epsmat=(NumericMatrix::diag(dd,1)*eps)/2;
#    NumericMatrix gmat(aa,dd);
#    for(int i=0;i<dd;i++){
#      for(int j=0;j<aa;j++){
#        NumericVector tempi=infcn(parvec + epsmat.column(i));
#        NumericVector tempii=infcn(parvec - epsmat.column(i));
#        gmat(j,i)=(tempi - tempii)[j]/eps;
#      }
#    }
#    return gmat;
#  }

## -----------------------------------------------------------------------------
# (g1,g2)=(x^3-y-1, x^3-y^2-1)
f1=function(x) c(x[1]^3-x[2]-1,x[1]^3-x[2]^2-1)
ini=c(1,2)
Gradmat(ini,f1)
Gradmat_cpp(ini,f1)
NRroot(ini,f1)

# (g1,g2,g3)=(x^4-y^2+z-1, x-y^2+z+3,x+y-3z-5)
f2=function(x) c(x[1]^4-x[2]^2+x[3]-1,x[1]-x[2]^2+x[3]+3,x[1]+x[2]-3*x[3]-5)
re=NRroot(c(1,2,3),f2)
print(re)

