#--------------------------------------------------------------- 
#-- delta.bspline2.R
#-- data generation function for calibration bias simulations 
#-- author: Sarah Hegarty
#-- date: 5 May 2025
#--------------------------------------------------------------- 


# B-spline Deviation Shifted
delta.curve <- function(z,delta = delta){
  basis <- splines:::bs(z, knots = c(0.15,0.35,0.55,0.97), degree = 3, Boundary.knots = range(z))
  
  
  coef <- as.matrix(c(delta,-delta,0,0,delta,0,0))
  basis %*% coef - delta
}

cal.curve <- function(z,delta = delta){
  cal.base(z) + delta.curve(z,delta)
}