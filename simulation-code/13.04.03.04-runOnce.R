# title: "13.04.03.04 - runOnce.R"
# author: SE Hegarty

runOnce <- function(funcs ="vc-funcs-2025-07-21.R"
                    , deltafunc = "delta.bspline1.R"
                    , n = c(1000,1000)
                    , sigma= c(0.5,0.5)
                    , delta = 0.5
                    , prefix = "13.04.03.04_temp-"
                    , M =1
                    , alpha = 0.05
                    , tau = c(0.5,0.9)
                    , tau.int= 0.1){
 
all.start <- Sys.time()

require(dplyr)
require(splines)
require(mgcv)
source(funcs)
source("cal.base.R")
source(deltafunc)
source("dataGen.R")


# ----------------------------------------------------------- #
# --------- Generate Data ----------------------------------- #
# ----------------------------------------------------------- #

df <- dataGen(n = n, sigma = sigma, delta= delta)

Y <- df$y.obs
s <- df$s
u <- df$risk
U <- unique(u)

# ----------------------------------------------------------- #
# --------- Bandwidth Selection ----------------------------- #
# ----------------------------------------------------------- #

# ---- Select pilot bandwidth for local cubic polynomial to minimize IRSC ---- #
rsc.start <- Sys.time()
# defining h grid using Fan and Gijbels approach:
U.max = max(U); U.min = min(U); n = length(U)
h.min <- (U.max - U.min)/n # this is too small for many u to invert all matrices so:
h.min <- max(h.min, 0.05*(U.max - U.min))
h.max <- (U.max - U.min)/2
h.grid <- seq(h.min, h.max, 0.02)
#h.grid <- c(0.3, 0.4, 0.5)
adj <- 0.7776

rsc.u0 <- NULL
irsc.h <- NULL
for(h in h.grid){
  rsc.h <- NULL
  for(x in U){
    r <- rsc(Y,s,U, u0=x, q = 1, a = 2, h=h)
    rsc.h <- rbind(rsc.h, c(x, h, r)) 
    colnames(rsc.h) <- c("u0","h","RSC")
  }
  irsc <- mean(rsc.h[,"RSC"])
  irsc.h <- rbind(irsc.h, c(h, irsc)) 
  colnames(irsc.h) <- c("h","IRSC")
  rsc.u0 <- rbind(rsc.u0, rsc.h)
}

h.star <- irsc.h %>% as.data.frame() %>% slice_min(IRSC, n = 1 ) %>% pull(h)
h.rsc <- h.star*adj

rsc.end <- Sys.time()
rsc.runtime <- rsc.end - rsc.start


# ---- Select local linear bandwidth by minimizing MSE(h) ---- #
bvs.start <- Sys.time()
# Get local cubic quantities necessary for MSE calculation
bvs <- list() 
for(j in 1:length(u)){
  bv.u <- zheng.mse.cubics(Y, s, u, u0 = u[j], h.star = h.rsc, q = 1, a = 2)
  bv.u
  bvs[[j]] <- bv.u
}
names(bvs) <- u
bvs.end <- Sys.time()
bvs.runtime <- bvs.end - bvs.start

mse.start <- Sys.time()
# Over grid of h, compute MSE
hgrid <- seq(0.02, 0.2, 0.02)
mse.h.u <- NULL
for(h in hgrid){
  for(u0 in u){
    m.u <- zheng.mse(Y = Y, s = s, u = u, u0 = u0, h = h, q=1, bvs = bvs)
    mse.h.u <- rbind(mse.h.u, c(h, u0, m.u[["mse"]], m.u[["mse.bias"]], m.u[["mse.var"]]))
    colnames(mse.h.u) <- c("h", "u0","MSE.u","MSE_bias.u","MSE_var.u")
  }
}

mse.h <- mse.h.u %>%
  as.data.frame() %>%
  group_by(h) %>%
  summarize(MSE.mean = mean(MSE.u)
            ,MSE.bias.mean = mean(MSE_bias.u)
            ,MSE.var.mean = mean(MSE_var.u))
mse.h.min <- mse.h %>%
  slice_min(MSE.mean, n = 1) 

h.mse <- mse.h.min %>% pull(h)

# Option 1: h.opt is h.mse
# h.opt <- h.mse

# Option 2: h.opt is 0.5*h.mse (approach used by Fan and Zhang (2000) in conjuction with setting bias = 0)
h.opt <- 0.5*h.mse 

# Option 3: (proposed) h.opt = argmin(MSE.bias) s.t. MSE <= MSE.min*1.10'
# min.MSE <- mse.h.min %>% pull(MSE.mean)
# h.opt <- mse.h %>% 
#             filter(MSE.mean <= min.MSE*1.10) %>% 
#             slice_min(MSE.bias.mean, n = 1) %>%
#             pull(h)

mse.end <- Sys.time()
mse.runtime <- mse.end - mse.start 


# ----------------------------------------------------------- #
# --------- Varying Coefficient Model ----------------------- #
# ----------------------------------------------------------- #

fit.start <- Sys.time()
# ---- Fit VC model with h.opt (linear) and h.rsc (cubic) ---- #
fit <- vcfit.bvs(Y,s,u, h=h.opt, h.bias = h.rsc, alpha = alpha)

ests <- fit %>%
  as.data.frame() %>%
  mutate(# calculate true curves
          true.int = cal.base(risk)
         ,true.s = delta.curve(risk, delta= delta)
         # Calculate D_a
         ,h.opt = h.opt
         ,D.int = D(alpha = alpha, h = h.opt, var = var.int)
         ,D.s = D(alpha = alpha, h = h.opt, var = var.s)
         # Calculate CIs for intercept
         ,lower.int = est.int - bias.int - D.int
         ,upper.int = est.int - bias.int + D.int
         # Calculate CIs for s
         ,lower.s = est.s - bias.s - D.s
         ,upper.s = est.s - bias.s + D.s
        )

# Global test statistic
l <- lh2(h.opt, c= 0, d = 1)
dnun <- d(h.opt, c= 0, d = 1)
sup <- ests %>%
  mutate(T.int = abs(est.int - bias.int)/sqrt(var.int)
         ,T.s = abs(est.s - bias.s)/sqrt(var.s)) %>%
  summarize(T.int = max(T.int)
            ,T.s = max(T.s)) 
global.T = l * (sup - dnun)
global.pT = 1 - FT(global.T)
 
# Global coverage
global.coverage <- ests %>%
  mutate(clower.int = if_else(lower.int <= true.int, 1, 0)
        ,cupper.int = if_else(upper.int >= true.int, 1, 0)
        ,cover.int = clower.int*cupper.int
        ,clowerband.int = min(clower.int)
        ,cupperband.int = min(cupper.int)
        ,coverage.int = min(cover.int)
        ,clower.s = if_else(lower.s <= true.s, 1, 0)
        ,cupper.s = if_else(upper.s >= true.s, 1, 0)
        ,cover.s = clower.s*cupper.s
        ,clowerband.s = min(clower.s)
        ,cupperband.s = min(cupper.s)
        ,coverage.s = min(cover.s)) %>%
      select(coverage.int, coverage.s) %>%
      unique()


# Local test statistic
local.T = local.pT = local.coverage = NULL
for(t in tau){
  t1 = tau - tau.int
  t2 = tau + tau.int
  ll <- lh2(h.opt, c= t1, d = t2)
  ldnun <- d(h.opt, c= t1, d = t2)
  lsup <- ests %>%
            filter(risk >= t1, risk <= t2) %>%
            mutate(T.int = abs(est.int - bias.int)/sqrt(var.int)
                   ,T.s = abs(est.s - bias.s)/sqrt(var.s)) %>%
            summarize(T.int = max(T.int)
                      ,T.s = max(T.s)) 
  lT = ll * (lsup - ldnun)
  local.T <- local.T %>% bind_rows(lT %>% as.data.frame() %>% mutate(tau = t))
  lpT = 1 - FT(local.T)
  local.pT <- local.pT %>% bind_rows(lpT %>% as.data.frame() %>% mutate(tau = t))
  
  # Local coverage
  lcoverage <- ests %>%
    filter(risk >= t1, risk <= t2) %>%
    mutate(clower.int = if_else(lower.int <= true.int, 1, 0)
           ,cupper.int = if_else(upper.int >= true.int, 1, 0)
           ,cover.int = clower.int*cupper.int
           ,clowerband.int = min(clower.int)
           ,cupperband.int = min(cupper.int)
           ,coverage.int = min(cover.int)
           ,clower.s = if_else(lower.s <= true.s, 1, 0)
           ,cupper.s = if_else(upper.s >= true.s, 1, 0)
           ,cover.s = clower.s*cupper.s
           ,clowerband.s = min(clower.s)
           ,cupperband.s = min(cupper.s)
           ,coverage.s = min(cover.s)) %>%
    select(coverage.int, coverage.s) %>%
    unique()
  local.coverage <- local.coverage %>% bind_rows(lcoverage %>% as.data.frame() %>% mutate(tau = t))
}

fit.end <- Sys.time()
fit.runtime <- difftime(fit.end, fit.start, units='mins')


vc.results <- list(ests = ests, h.mse = h.mse, h.rsc = h.rsc, h.opt = h.opt
                   ,global.T = global.T, global.pvalues = global.pT, global.coverage= global.coverage
                   ,local.taus = tau, local.interval = tau.int
                   ,local.T= local.T, local.pvalues = local.pT, local.coverage = local.coverage
                   , rsc.runtime = rsc.runtime
                   , bvs.runtime = bvs.runtime, mse.runtime = mse.runtime
                   ,vc.runtime = fit.runtime
                   ,mse.h = mse.h, irsc.h = irsc.h)

# ----------------------------------------------------------- #
# --------- Comparative Methods ----------------------------- #
# ----------------------------------------------------------- #

# ------------- LRT  -----------------------------------------#
lrt.start <- Sys.time()
lm.fit.full <- lm(y.obs ~ risk*s, data = df)
lm.fit.red <- lm(y.obs ~ risk, data = df)

# based on chi-square test
Tlm <- anova(lm.fit.red, lm.fit.full, test = "LRT")
lrt.stop <- Sys.time()
lrt.runtime <- difftime(lrt.stop, lrt.start, units='secs')
lrt.results <- list(T.s = Tlm[2,4], df = Tlm[2,3], pvalue = Tlm[2,5], lrt.runtime = lrt.runtime)


# ------------- Additive Model: spline + linear ------------- #
aml.start <- Sys.time()
aml.fit <- with(df, mgcv:::gam(y.obs ~ s(risk, bs= 'cr') + s + s:risk))
aml.sum <- summary(aml.fit)
aml.cov <- aml.fit[["Ve"]]
aml.beta <- aml.fit[["coefficients"]]
rownames(aml.cov) = colnames(aml.cov) <- names(aml.beta)
# sqrt(diag(aml.cov)[1:3]) # matches std.error of parameters estimates in summary(aml.fit)
lincomps <- c("s","s:risk")
Vinv <- solve(aml.cov[lincomps, lincomps])
Taml <- as.numeric( crossprod(aml.beta[lincomps], Vinv) %*% aml.beta[lincomps] )
p.aml <- 1 - pf(Taml, df1 = length(lincomps), df2 = aml.sum[["n"]] - aml.sum[["np"]])
# aml.lincomp <- function(x){
#           aml.beta[2] + aml.beta[3]*x
# }
# curve(aml.lincomp, ylim = c(-0.1,0.1))  
# abline(h=0, lty=2)
aml.stop <- Sys.time()
aml.runtime <- difftime(aml.stop, aml.start, units='secs')
aml.results <-list(T.s = Taml, ndf = length(lincomps), ddf=  aml.sum[["n"]] - aml.sum[["np"]], 
                   pvalue =p.aml, aml.fit = aml.fit, aml.runtime = aml.runtime)

# ------------- Additive Model: spline + spline ------------- #
ams.start <- Sys.time()
ams.fit <- with(df, mgcv:::gam(y.obs ~ s(risk, bs= 'cr') + s(delta, bs='cr')))
ams.results <- anova(ams.fit)
Tams <- as.data.frame(ams.results$s.table)
ams.stop <- Sys.time()
ams.runtime <- difftime(ams.stop, ams.start, units='secs')
ams.results <-list(T.s = Tams[2,3], pvalue = Tams[2,4] ,T.full = Tams,  ams.runtime = ams.runtime)

# ------------- Threshold Calibration Error ------------- #
TCE.fit = vcfit.bvs(Y, s, u, h = h.opt, h.bias = h.rsc, alpha= 0.05, ez = tau)

all.end <- Sys.time()
all.runtime <-difftime(all.end, all.start, units='mins')


results <- list(vc = vc.results
          , lrt = lrt.results
          , aml = aml.results
          , ams = ams.results
          , TCE = TCE.fit
          , all.runtime = all.runtime)

save(results, file = paste0("results/",prefix,M,".RData"))


return(results)
}
