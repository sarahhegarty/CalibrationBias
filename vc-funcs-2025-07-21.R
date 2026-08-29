#--------------------------------------------------- 
#-- vc-funcs-2025-07-21.R
#-- support functions for varying coefficient model
#-- implementation to estimate calibration bias
#-- author: Sarah Hegarty
#-- date: 21 July 2025
#--------------------------------------------------- 


colMax <- function(data) sapply(data, max, na.rm = TRUE)

# .kernel is copied from tvReg package
.kernel <- function(x, bw, tkernel = "Triweight") {
  x <- x/c(bw)
  value <- numeric(length(x))
  if (tkernel == "Gaussian") {
    value <- exp(-0.5 * x^2)/2.506628274631
  }
  else if (tkernel == "Epa") {
    index <- (x >= -1 & x <= 1)
    value[index] <- 0.75 * (1 - x[index]^2)
  }
  else if (tkernel == "Triweight") {
    index <- (x >= -1 & x <= 1)
    value[index] <- 1.09375 * (1 - x[index]^2)^3
  }
  return(value)
}

FT <- function(t){
  exp(-2*exp(-t))
}

lh2 <- function(h,c = 0, d =1){
  lh = (-2 * log(h/(d-c)))**(1/2)
  return(lh)
}

d <- function(h,c=0,d=1){
  nu0 = 0.6
  K2= 1.5
  lh = lh2(h,c,d)
  d= lh + 1/lh * log( 1/(4*nu0*pi) * K2 )
  return(d)
}

D <- function(alpha, h, var, c=0, d=1){
  # from p 722 of Fan and Zhang 2000
 
  lh = lh2(h,c,d)
  
  d = d(h,c,d)
  clogalpha = log(-log(1-alpha))
  mult = d + ( log(2) - clogalpha ) / lh
  
  D = mult * sqrt(var)
  return(D)
}


rsc <- function(Y, s, u, u0, q = 1, a = 2, h){
  
  df = data.frame(Y,s,u)
  #ez = unique(u)
  
  # Setup
  x = cbind(1,s)
  U = u
  Uu = U - u0
  
  p <- ncol(x)
  
  
  # ---------------- Data Setup ------------------ #
  # with kernel matrix with current bandwidth
  K = .kernel(Uu, bw = h , tkernel = "Epa") 
  K.index = which(K != 0)
  W = diag(K[K.index])
  
  # construct design matrix for local polynomial of order q+a
  xtemp = x[K.index, ]
  ytemp = Y[K.index]
  
  X = NULL 
  for(k in 1:p){ # loop over number of variables
    for(j in 0:(q+a)){ # loop over order of polynomial
      X = cbind(X, xtemp[,k] * (Uu[K.index]^j) )
    }
  }
  
  S = crossprod(X,W) %*% X
  S.star = crossprod(X,W**2) %*% X
  S.inv = solve(S)
  V = (S.inv %*% S.star %*% S.inv)[1,1]
  
  # Get kernel weighted polynomial parameter estimates
  bhat <- S.inv %*% t(X) %*% W %*% ytemp
  
  # Get estimated y
  yhat <- X %*% bhat
  
  # Get sigma hat
  trace <- sum(diag(W)) - sum(diag(S.inv %*% S.star))
  sigma.hat <- 1/trace * sum( (ytemp - yhat)**2 * diag(W))
  
  # Calculate RCS
  rsc <- sigma.hat * (1 + (q + a + 1)*V)
  return(rsc)
}

zheng.mse.cubics <- function(Y,s,u,u0,h.star, q= 1,a = 2){
  df = data.frame(Y,s,u)
  #ez = unique(u)
  
  # Setup
  x = cbind(1,s)
  U = u
  Uu = U - u0
  
  p <- ncol(x)
  
  # Local polynomial of order q + a for bias, variance, covariance estimates
  K.star = .kernel(Uu, bw = h.star , tkernel = "Epa") 
  K.star.index = which(K.star != 0)
  W.star = diag(K.star[K.star.index])
  
  # construct design matrix for local polynomial of order q+a
  x.star = x[K.star.index, ]
  y.star = Y[K.star.index]
  
  X.star.full = NULL 
  for(k in 1:p){ # loop over number of variables
    for(j in 0:(q+a)){ # loop over order of polynomial
      X.star.full = cbind(X.star.full, x[,k] * (Uu^j) )
    }
  }
  
  X.star <- X.star.full[K.star.index, ]
  S.star = crossprod(X.star,W.star) %*% X.star
  S2.star = crossprod(X.star,W.star**2) %*% X.star
  S.inv.star = solve(S.star)
  VV.star = (S.inv.star %*% S2.star %*% S.inv.star)
  
  # Get kernel weighted polynomial parameter estimates
  bhat.star <- S.inv.star %*% t(X.star) %*% W.star %*% y.star
  
  # Get estimated y
  yhat.star <- X.star %*% bhat.star
  
  # Variance Estimate (sigma.hat)
  ESS <- sum(diag(W.star)) - sum(diag(S.inv.star %*% S2.star))
  sigma.hat <- 1/ESS * sum( (y.star - yhat.star)**2 * diag(W.star))
  
  # Covariance Calculation (Omega)
  V.star = NULL 
  for(j in 0:(q+a)){ # loop over order of polynomial
    V.star = cbind(V.star, Uu[K.star.index]^j)
  }
  
  VW <- crossprod(V.star, W.star)
  VWV <- solve(VW %*% V.star) %*% VW
  
  Z11 <- x.star[,1]*x.star[,1]
  Z12 <- x.star[,1]*x.star[,2]
  Z21 <- x.star[,2]*x.star[,1]
  Z22 <- x.star[,2]*x.star[,2]
  e1g <- rep(0,q+a+1); e1g[1] <-1
  r11 <- e1g %*% VWV %*% Z11
  r12 <- e1g %*% VWV %*% Z12
  r21 <- e1g %*% VWV %*% Z21
  r22 <- e1g %*% VWV %*% Z12
  Omega <- matrix(c(r11, r12, r21, r22), nrow = 2, ncol=2, byrow = TRUE)
  
  # Bias calculation (Xd)
  bias.terms <- c(3,4,7,8)
  dhat <- bhat.star[bias.terms]
  Xd <-  X.star.full[,bias.terms] %*% dhat
  
  cubics <- list(Xd = Xd, Omega = Omega, sigma.hat = sigma.hat)
  return(cubics)
}

zheng.mse <- function(Y, s, u, u0, h, q = 1, bvs){
  # browser()
  df = data.frame(Y,s,u)
  #ez = unique(u)
  
  # Setup
  x = cbind(1,s)
  U = u
  Uu = U - u0
  
  p <- ncol(x)
  
  # ---------------- Data Setup ------------------ #
  # with kernel matrix with current bandwidth
  K = .kernel(Uu, bw = h , tkernel = "Epa") 
  K.index = which(K != 0)
  W = diag(K[K.index])
  
  # construct design matrix for local polynomial of order q+a
  xtemp = x[K.index, ]
  ytemp = Y[K.index]
  
  X = NULL 
  for(k in 1:p){ # loop over number of variables
    for(j in 0:q){ # loop over order of polynomial
      X = cbind(X, xtemp[,k] * (Uu[K.index]^j) )
    }
  }
  
  S = crossprod(X,W) %*% X
  S2 = crossprod(X,W**2) %*% X
  S.inv = solve(S)
  VV = (S.inv %*% S2 %*% S.inv)
  
  # Get kernel weighted polynomial parameter estimates
  bhat.q <- S.inv %*% t(X) %*% W %*% ytemp
  
  # Get estimated y from linear
  yhat.q <- X %*% bhat.q
  
  # Calculate MSE
  bv.u <- bvs[[paste0(u0)]]
  Xd <- bv.u[["Xd"]][K.index, ]
  Omega <- bv.u[["Omega"]]
  sigma.hat <- bv.u[["sigma.hat"]]
  e1q <- e.vec(1,2)
  EE <- tcrossprod(e1q,e1q) 
  mse1 <- crossprod(Xd, W) %*% X %*% S.inv %*% (Omega %x% EE)%*% S.inv %*% t(X) %*% W %*% Xd
  mse2 <- sum(diag(VV %*% (Omega %x% EE))) * sigma.hat
  mse <-as.numeric( mse1 + mse2)
  mse.out <- list(mse = mse, mse.bias = mse1, mse.var = mse2)
  return(mse.out)
} 

e.vec <- function(j,q){
  e.q <- numeric(q)
  e.q[j] <- 1
  return(e.q)
}

vcfit.bvs <- function(Y, s, u, h , h.bias, alpha =0.05, ez = NULL){
  #browser()
  df = data.frame(Y,s,u)
  if(is.null(ez)){ez = unique(u)}
  
  require(dplyr)
 
  n = nrow(df)
  
  # ---------------- Asymptotic-based estimators  -------------------- #
  # Setup
  x = cbind(1,s)
  U = u
  
  p = ncol(x)
  n = nrow(x)
  q = 1 # local linear polynomial
  a = 2 # q+a= local cubic polynomial
  
  # utilities for local linear fit
  kappa = p*(q+1)
  mainfx <- seq(1,kappa,q+1)
  
  # initialize output vector
  vcest <- matrix(0,nrow = length(ez), ncol = 3 * p + 2)
  colnames(vcest) <- c("risk"
                       ,"est.int","est.s"
                       ,"bias.int","bias.s"
                       ,"var.int","var.s"
                       ,"sigma2.hat"
  )
  
  for(r in 1:length(ez)){
    u0 = ez[r]
    Uu = U - u0
    # ---------------- Data Setup ------------------ #
    # construct kernel matrices with local linear bandwidth
    K = .kernel(Uu, bw = h , tkernel = "Epa") 
    K.index = which(K != 0)
    W = diag(K[K.index])
    Uu = Uu[K.index]
    
    # construct design matrix for local polynomial of order q
    xtemp = x[K.index, ]
    
    X = NULL 
    for(k in 1:p){ # loop over number of variables
      for(j in 0:q){ # loop over order of polynomial
        X = cbind(X, xtemp[,k] * (Uu**j) )
      }
    }
    
    # ------ Weighted Least Squares Estimate of Local Order "q" Polynomial with h.opt bw------ #
    XWXinv <- solve((crossprod(X,W) %*% X))
    H <- XWXinv %*% crossprod(X,W) 
    beta <- H %*% Y[K.index]
    ahat.u <- beta[mainfx]
    
    # ------ Get BVS Estimates from "q + a" Order Polynomial with hRSC bandwidth ------ #
    bvs <- zheng.mse.cubics(Y, s, u, u0 = u0, h.star = h.bias, q = q, a = a)
    tau.u <- bvs[["Xd"]][K.index,]
    sigma.hat <- bvs[["sigma.hat"]]
    
    # estimate covariance per Fan and Zhang 2008, p 182 middle
    XW2X = crossprod(X,W**2) %*% X
    covhat = XWXinv %*% XW2X %*%XWXinv
    e.1 <- e.vec(1,4)
    e.3 <- e.vec(3,4)
    IO <- rbind(e.1,e.3)
    cov.hat <- IO %*% XWXinv %*% XW2X %*%XWXinv %*% t(IO) * sigma.hat
    var.hat <- diag(cov.hat)
    
    # estimate bias per Fan and Zhang 2008, p 182 top
    bias.hat = IO %*% XWXinv %*% t(X) %*% W %*% tau.u 
    
    # output results
    vcest[r,] <- c(u0, ahat.u, bias.hat, var.hat, sigma.hat)
    
  }
  return(vcest)
}
    

