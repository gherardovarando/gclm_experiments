args = commandArgs(trailingOnly=TRUE)  


lower <- FALSE
p <- 10
Ps <- c(10)
rep <- 100
nlambda <- 100
N <- 1000
bpath <- "simulations/"
if (length(args) != 0){ 
  message("arguments found, parsing...")
  la <- length(args) 
  ixpath <- which(args %in% "path")
  if (length(ixpath) == 1){
     if (ixpath < la){
         bpath <- args[ixpath + 1] 
         message("base path set to ", bpath)
     }
  }
  ixp    <- which(args %in% "p")  
  if (length(ixp) == 1){
     if (ixp < la){
         p <- as.numeric(args[ixp + 1]) 
         message("p set to ", p)
     }
  }
  ixN <- which(args %in% "N")
  if (length(ixN) == 1){
     if (ixp < la){
         N <- as.numeric(args[ixN + 1]) 
         message("N set to ", N)
     }
  }
  ixP    <- which(args %in% "P") 
  if (length(ixP) == 1){
     if (ixP < la){
         if (ixP > max(ixp, ixpath)){
             Ps <- as.numeric(args[(ixP + 1):(la)]) 
         }else{
             Ps <- as.numeric(args[ixP + 1])
         }
         message("P(s) set to ", Ps)
     }
  }
}

library(gclm)
library(igraph)
library(ggplot2)
source("functions/util.R")

message("packages and util functions loaded correctly")

#lambdaseq <- exp( (1 - (100:1)) / 10)
lambdaseq <- 10^(seq(-4,0,length = 100))
message("starting experiments..")
for (P in Ps) {
  for (k in c(1,2,3,4)) {
    path <- paste0(bpath,
                   "/p",
                   p,
                   "/P",
                   P,
                   "/k",
                   k,
                   "/N",
                    N,
                    "/")
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    d <- k / p
    for (r in 1:rep) {
      Btrue <-
        rStableMetzler(
          n = P,
          p = d,
          lower = lower,
          rfun = function(n)
            rnorm(n, sd = 1),
          rdiag = rnorm
        )
      Ctrue <- diag(runif(P))
      exper <- rOUinv(N, Btrue, C = Ctrue)
      results <- list()
      times <- list()
      C0 <- diag(p)
      Sigmahat <- cor(exper$data[1:N,1:p])
      Bstart <- -0.5 *  solve(Sigmahat)
      message("running loglik-inf")
        tllb <-
          system.time(
            resllb <- gclm.path(
              Sigmahat,
              eps = 1e-4,
              B = Bstart,
              maxIter = 100,
              job = 0,
              lambdas =  6*lambdaseq,
              lambdac = -1
            )
          )

      message("running frobenius-inf")
        tfrobenius <-
          system.time(
            resfrobenius <- gclm.path(
              Sigmahat,
              eps = 1e-4,
              B = Bstart,
              maxIter = 100,
              job = 0,
              lambdas = 6*lambdaseq,
              lambdac = -1,
              loss = "frobenius"
            )
          )
      message("running loglik-0.01")
        tpnll <-
          system.time(
            respnll <- gclm.path(
              Sigmahat,
              eps = 1e-4,
              B = Bstart,
              maxIter = 100,
              job = 0,
              lambdas = 6*lambdaseq,
              lambdac = 0.01
            )
          )
        message("runnign lasso")
        tlasso <-
          system.time(reslasso <- lassoB(Sigmahat, C = C0,
                                         lambda = lambdaseq))

        message("runnign glasso")
        tglasso <- system.time(resglasso <- glassoB(Sigmahat,
                                                    lambda = lambdaseq * 
                                                    max(diag(Sigmahat))))
        message("runnign covthr")
        tcovthr <- system.time(rescovthr <- covthr(Sigmahat))
        times[[paste0(N)]] <-
          list(loglik = tllb,
               frobenius = tfrobenius,
               pnll = tpnll,
               lasso = tlasso,
               glasso = tglasso,
               covthr = tcovthr)
        results[[paste0(N)]] <- list(loglik = resllb,
                                     frobenius = resfrobenius,
                                     pnll = respnll,
                                     lasso = reslasso,
                                     glasso = resglasso,
                                     covthr = rescovthr)
      name <- paste0("rep", r, ".RData")
      message("DONE rep ", r, "P = ", P,
              " p=", p, " k=", k)
      save(
        file = paste0(path, name),
        list = c(
          "results",
          "exper",
          "times",
          "p",
          "d",
          "P",
          "N",
          "lower"
        )
      )
    }
  }
}
