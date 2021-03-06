## ----------------------------------------------------------------------------
## Name : trapint.r
## ----------------------------------------------------------------------------
## Authors: Matt P. Wand
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description: Trapezoidal integration
## ----------------------------------------------------------------------------

trapint <- function(xgrid,fgrid) {
   ng <- length(xgrid)
   xvec <- xgrid[2:ng] - xgrid[1:(ng-1)]
   fvec <- fgrid[1:(ng-1)] + fgrid[2:ng]
   integ <- sum(xvec*fvec)/2
   return(integ)
}
