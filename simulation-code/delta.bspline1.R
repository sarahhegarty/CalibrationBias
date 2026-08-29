#--------------------------------------------------------------- 
#-- delta.bspline1.R
#-- data generation function for calibration bias simulations 
#-- author: Sarah Hegarty
#-- date: 15 August 2025
#--------------------------------------------------------------- 

# Define B-spline for difference
delta.curve <- function(z, delta = delta){
  require(splines)
  basis <- bs(z, knots = c(0.15,0.35,0.55,0.97), degree = 3, Boundary.knots = c(0,1))
  
  coef <- as.matrix(c(0,delta,-delta,-delta,delta,0,0))
  basis %*% coef
}

# Get calibration curve for group with deviation
cal.curve <- function(z,delta =delta){
  cal.base(z) + delta.curve(z,delta )
}

