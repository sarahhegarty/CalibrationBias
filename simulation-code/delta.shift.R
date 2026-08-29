#--------------------------------------------------------------- 
#-- delta.shift.R
#-- data generation function for calibration bias simulations 
#-- author: Sarah Hegarty
#-- date: 5 May 2025
#--------------------------------------------------------------- 

delta.curve <- function(z, delta = delta){
  delta
}


cal.curve <- function(z,delta = delta){
  cal.base(z) + delta.curve(z,delta)
}