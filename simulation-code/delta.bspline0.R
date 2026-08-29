#--------------------------------------------------------------- 
#-- delta.bspline0.R
#-- data generation function for calibration bias simulations 
#-- author: Sarah Hegarty
#-- date: 31 July 2025
#--------------------------------------------------------------- 

# Define B-spline for difference
delta.curve <- function(z, delta = delta){
  basis <- bs(z, knots = c(0.15,0.35,0.55,0.97), degree = 3, Boundary.knots = range(z))
  
  coef <- as.matrix(c(-delta,-delta,0,0,delta,0.5*delta,0.5*delta))
  basis %*% coef
}

# Get calibration curve for group with deviation
cal.curve <- function(z,delta =delta){
  cal.base(z) + delta.curve(z,delta )
}

