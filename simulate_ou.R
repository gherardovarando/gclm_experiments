library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")


rescaleC <- FALSE
lower <- FALSE
p <- 10
Ps <- c(10, 12, 15, 20)
rep <- 100
nlambda <- 100
lambdaseq <- c(0, exp(10 * (1 - (nlambda:1)) / nlambda))
for (P in Ps) {
  for (k in c(1,2,3,4)) {
    path <- paste0("simulations/",
                   "/p",
                   p,
                   "P",
                   P,
                   "/k",
                   k,
                   "/")
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    d <- k / P
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
      Ctrue <- diag(runif(p))
      exper <- rOUinv(5000, Btrue, C = Ctrue)
      results <- list()
      times <- list()
      for (N in c(50,
                  100,
                  200,
                  300,
                  400,
                  500,
                  600,
                  700,
                  800,
                  900,
                  1000,
                  2000,
                  3000,
                  4000,
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
              eps = 1e-6,
              C = C0,
              maxIter = 5000,
              job = 11,
              lambdas = 3 * lambdaseq
            )
          )
        tlasso <-
          system.time(reslasso <- lassoB(Sigmahat, C = C0,
                                         lambda = lambdaseq))
        tglasso <- system.time(resglasso <- glassoB(Sigmahat,
                                                    lambda = lambdaseq * max(diag(Sigmahat))))
        times[[paste0(N)]] <-
          list(loglik = tllb,
               lasso = tlasso,
               glasso = tglasso)
        results[[paste0(N)]] <- list(loglik = resllb,
                                     lasso = reslasso,
                                     glasso = resglasso)
      }
      name <- paste0("rep", r, ".RData")
      message("DONE rep ", r, typeC, typeB,
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
