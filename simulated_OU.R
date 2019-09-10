library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")


p <- 10
rep <- 100
lower <- FALSE

for (k in c(2)){
  d <- k / p
  for (N in c(50, 100, 200, 300, 400, 500, 1000)){
    dir.create("simulated_OU", showWarnings = FALSE)
    path <- paste0("simulated_OU/identityC/p", p, "/k", k, "/N", N)
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    sapply(1:rep, function(i){
      Btrue <- rStableMetzler(n = p, p = d, lower = lower, rfun = function(n) 
        rnorm(n, sd = 1), 
        rdiag = rnorm)
      #Ctrue <- diag(runif(p))
      Ctrue <- diag(p)
      exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))
#      Sigmahat <- cov(exper$data)
      Sigmahat <- cor(exper$data)
      resllb <- llBpath(Sigmahat, eps = 1e-12, 
                         maxIter = 5000, job = 11, 
                         lambdas = seq(0,3,length.out = 100))
      reslsb <- lsBpath(Sigmahat, eps = 1e-15, 
                        lambdas = seq(0,3,length.out = 100),
                        maxIter = 5000, job =11)
      reslasso <- lassoB(Sigmahat)
      results <- list(loglikL1B = resllb,
                      lsL1B = reslsb,
                      lasso = reslasso)
      name <- paste0("rep", i,".RData")
      save(file = paste0(path, "/", name), list = c("results",
                                                    "exper", 
                                                    "p", 
                                                    "N", "d", "lower"))
    })  
  }
}

