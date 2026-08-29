#--------------------------------------------------------------- 
#-- dataGen.R
#-- data generation function for calibration bias simulations 
#-- author: Sarah Hegarty
#-- date: 5 May 2025
#--------------------------------------------------------------- 

dataGen <- function(n = c(100,100)
                         ,sigma= c(0.5,0.5)
                          ,delta=0){
  
  # group 0 - with base curve
  s = rep(0,n[1])
  risk = runif(n[1])
  y.pred = cal.base(risk)
  epsilon = rnorm(n[1],0,sigma[1])
  y.obs = y.pred + epsilon
  df1 <- data.frame(s, risk, y.pred, epsilon, y.obs)
  
  # group 1 - with delta curve
  s = rep(1,n[2])
  risk = runif(n[2])
  y.pred = cal.curve(risk, delta=delta)
  epsilon = rnorm(n[2],0,sigma[2])
  y.obs = y.pred + epsilon
  df2 <- data.frame(s, risk, y.pred, epsilon, y.obs)
  
  df <- df1 %>% bind_rows(df2) %>%
    mutate(delta = risk*s)
  
  return(df)
}
