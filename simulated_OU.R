library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")


p <- 10
rep <- 100
lower <- FALSE
for (k in c(1, 2, 3, 4)){
  path <- paste0("simulated_OU/randomC/p", p, "/k", k, "/")
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  d <- k / p
  for (r in 1:rep){
  Btrue <- rStableMetzler(n = p, p = d, lower = lower, rfun = function(n) 
    rnorm(n, sd = 1), 
    rdiag = rnorm)
  Ctrue <- diag(runif(p))
  exper <- rOUinv(5000, Btrue, C = Ctrue)
  results <- list()
  times <- list()
  for (N in c(50, 100, 200, 300, 400, 500,
              600, 700, 800, 900, 
              1000, 2000, 3000, 4000, 5000)){
    Sigmahat <- cor(exper$data[1:N,])
    tllb <- system.time(resllb <- llBpath(Sigmahat, eps = 1e-6, C = diag(p), 
                      maxIter = 5000, job = 11, 
                      lambdas = seq(0,3,length.out = 100)))
    tlsb <- system.time(reslsb <- lsBpath(Sigmahat, eps = 1e-6, C = diag(p),
                      lambdas = seq(0,3,length.out = 100),
                      maxIter = 5000, job =11))
    tlasso <- system.time(reslasso <- lassoB(Sigmahat, C = diag(p), 
                       lambda = seq(0.01, 0, length.out = 100)))
    tglasso <- system.time(resglasso <- glassoB(Sigmahat, 
                         lambda = seq(max(abs(Sigmahat)) / 20,
                                      max(abs(Sigmahat)), length.out = 100)))
    times[[paste0(N)]] <- list(loglik = tllb, frobenius = tlsb, lasso = tlasso, 
                               glasso = tglasso)
    results[[paste0(N)]] <- list(loglik = resllb,
                    frobenius = reslsb,
                    lasso = reslasso,
                    glasso = resglasso)
  }
  name <- paste0("rep", r,".RData")
  save(file = paste0(path, name), list = c("results",
                                                "exper",
                                           "times",
                                                "p", 
                                                "d", "lower"))
}
}

