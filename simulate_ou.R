args = commandArgs(trailingOnly=TRUE)  

library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

message("packages and util functions loaded correctly")

rescaleC <- FALSE
lower <- FALSE
p <- 10
Ps <- c(10, 12, 15, 20, 25, 30, 35, 40)
rep <- 100
nlambda <- 100
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
lambdaseq <- c(exp(10 * (1 - (nlambda:1)) / nlambda))
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
      exper <- rOUinv(5000, Btrue, C = Ctrue)
      results <- list()
      times <- list()
      for (N in c(
                 # 50,
                  100,
                #  200,
                #  300,
                #  400,
                  500,
               #   600,
               #   700,
               #   800,
               #   900,
                  1000,
               #   2000,
               #   3000,
               #   4000,
                  5000)) {
        if (rescaleC) {
          C0 <- diag(p)
          Sigmahat <- cov(exper$data[1:N,1:p])
          diag(C0) <- diag(C0) / diag(Sigmahat)
          Sigmahat <- cov2cor(Sigmahat)
        } else{
          C0 <- diag(p)
          Sigmahat <- cor(exper$data[1:N,1:p])
        }
        tllb <-
          system.time(
            resllb <- llBpath(
              Sigmahat,
              eps = 1e-4,
              C = C0,
              maxIter = 100,
              job = 11,
              lambdas = 3 * lambdaseq
            )
          )
#        tfrobenius <-
#          system.time(
#            resfrobenius <- lsBpath(
#              Sigmahat,
#              eps = 1e-4,
#              C = C0,
#              maxIter = 100,
#              job = 11,
#              lambdas = 3 * lambdaseq
#            )
#          )
        tlasso <-
          system.time(reslasso <- lassoB(Sigmahat, C = C0,
                                         lambda = lambdaseq))
#        tlassoc <-
#          system.time(reslassoc <- lassoB(cov(exper$data[1:N,1:p]), 
#                                          C = Ctrue[1:p,1:p],
#                                         lambda = lambdaseq))
        tglasso <- system.time(resglasso <- glassoB(Sigmahat,
                                                    lambda = lambdaseq * 
                                                    max(diag(Sigmahat))))
        tcovthr <- system.time(rescovthr <- covthr(Sigmahat))
        times[[paste0(N)]] <-
          list(loglik = tllb,
             #  frobenius = tfrobenius,
               lasso = tlasso,
             #  lassoc = tlassoc,
               glasso = tglasso,
               covthr = tcovthr)
        results[[paste0(N)]] <- list(loglik = resllb,
                                   #   frobenius = resfrobenius,
                                     lasso = reslasso,
                                   #  lassoc = reslassoc,
                                     glasso = resglasso,
                                     covthr = rescovthr)
      }
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
          "lower",
          "rescaleC"
        )
      )
    }
  }
}
