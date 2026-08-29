library(dplyr)
library(microbenchmark) 
library(parallel) 
library(foreach)
library(iterators) 
library(doParallel)

# --- inputs --- #
prefix = "13.04.03.04_bs0-0.2-n500-"
M = 1500
seed = 202508152
cores = 20

  set.seed(seed)
  
  seeds <- sample(1:999999999, M, replace = FALSE)
  
  # setup parallelization
  maxCores <- detectCores()
  maxCores
  ncluster <- min(maxCores-1,cores)
  simCluster <- makeCluster(ncluster, type = "PSOCK")
  registerDoParallel(simCluster)
  mcStart <- Sys.time()
  
  # parallelize over independent MC replicates
  sum.M <- foreach(m = 1:M, .inorder = FALSE, .combine=rbind
                   , .packages=c("dplyr","tidyr","mgcv","splines") ) %dopar% {
                     
            source("13.04.03.04-runOnce.R") 
            set.seed(seeds[m])          
            r <- runOnce(funcs ="vc-funcs-2025-07-21.R"
                         , deltafunc = "delta.bspline0.R"
                         , n = c(500,500)
                         , sigma=c(0.5,0.5)
                         , delta = 0.2
                         , prefix = "13.04.03.04_bs0-0.2-n500-"
                         , M = m
                         , alpha = 0.05 
                         , tau = c(0.5,0.9)
                         , tau.int = 0.1)
            
            
           # extract primary results #
            vc <- r[["vc"]]
            lrt <- r[["lrt"]]
            aml <- r[["aml"]]
            ams <- r[["ams"]]
            TCE <- r[["TCE"]]
            rt <- r[["all.runtime"]]
            
            out <-list(m=m
                      , vc.p_global = vc[["global.pvalues"]][["T.s"]]
                      , vc.cover_global = vc[["global.coverage"]][["coverage.s"]]
                      , vc.runtime= vc[["vc.runtime"]]
                      , vc.local_t1 = vc[["local.interval"]][1]
                      , vc.local_t2 = vc[["local.interval"]][2]
                      , vc.p_local = vc[["local.pvalues"]][["T.s"]]
                      , vc.cover_local = vc[["local.coverage"]][["coverage.s"]]
                      , lrt.p = lrt[["pvalue"]]
                      , lrt.runtme = lrt[["lrt.runtime"]]
                      , aml.p = aml[["pvalue"]]
                      , aml.runtime = aml[["aml.runtime"]]
                      , ams.p = ams[["pvalue"]]
                      , ams.runtime = ams[["ams.runtime"]]
                      , total.runtime = rt[[1]]
                      )
          
            outvec <- unlist(out)
            return(outvec)
        }
  # end foreach
  
  # close out parallelization
  stopCluster(simCluster)
  mcStop <- Sys.time()
  mcTime <- mcStop-mcStart
  
  save(sum.M,file=paste0(prefix,"MC_",M,".Rdata"))
  
  # Run time 
  mcTime
  
