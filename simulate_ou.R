library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")


p <- 10
rep <- 100
nlambda <- 100
for (randomC in c(TRUE, FALSE)){
    for (lower in c(FALSE, TRUE)){
      typeC <- ifelse(randomC, "randomC", "identityC")
      typeB <- ifelse(lower, "lowertriangular", "general")
      for (k in c(1, 2, 3, 4)){
        path <- paste0("simulations/", typeC, "/", typeB, "/p", p, "/k", k, "/")
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
        d <- k / p
        for (r in 1:rep){
          Btrue <- rStableMetzler(n = p, p = d, lower = lower, rfun = function(n) 
            rnorm(n, sd = 1), 
            rdiag = rnorm)
          if (randomC){
            Ctrue <- diag(runif(p))
          }else{
            Ctrue <- diag(p) 
          }
          C0 <- diag(p)
          exper <- rOUinv(5000, Btrue, C = Ctrue)
          results <- list()
          times <- list()
          for (N in c(50, 100, 200, 300, 400, 500,
                      600, 700, 800, 900, 
                      1000, 2000, 3000, 4000, 5000)){
            Sigmahat <- cov(exper$data[1:N,]) * (N - 1) / N
            #D <- diag( 1 / sqrt(diag(Sigmahat)))
            #C0 <- D %*% C0 %*% D
            Sigmahat <- cov2cor(Sigmahat)
            tllb <- system.time(resllb <- llBpath(Sigmahat, eps = 1e-6, C = C0, 
                                                  maxIter = 5000, job = 11, 
                                                  lambdas = seq(0,3,length.out = nlambda)))
            tlsb <- system.time(reslsb <- lsBpath(Sigmahat, eps = 1e-6, C = C0,
                                                  lambdas = seq(0,3,length.out = nlambda),
                                                  maxIter = 5000, job = 11))
            tlasso <- system.time(reslasso <- lassoB(Sigmahat, C = C0, nlambda = nlambda))
            tglasso <- system.time(resglasso <- glassoB(Sigmahat, 
                                                        lambda = seq(0,
                                                                     max(abs(Sigmahat)), length.out = nlambda)))
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
    }
}
