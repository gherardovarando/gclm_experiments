library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")


p <- 10
rep <- 100

for (k in c(2)){
  d <- k / p
  for (N in c(50, 100, 200, 300, 400, 500, 1000)){
    dir.create("simulated_OU", showWarnings = FALSE)
    dirname <- paste0("p",p,"N",N)
    dir.create(paste0("simulated_OU/",dirname), showWarnings = FALSE)
    path <- paste0("simulated_OU/", dirname, "/", "d", d, "/")
    dir.create(path, showWarnings = FALSE)
    sapply(1:rep, function(i){
      Btrue <- rStableMetzler(n = p, p = d, lower = FALSE, rfun = function(n) 
        rnorm(n, sd = 1), 
        rdiag = rnorm)
      Ctrue <- diag(runif(p))
      exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))
      Sigmahat <- cor(exper$data)
      results <- llBpath(Sigmahat, eps = 1e-12, 
                         maxIter = 5000, job = 11, 
                         lambdas = seq(0,4,length.out = 100))
      roc <- t(sapply(results, function(res){
        ix <- lower.tri(Btrue) | upper.tri(Btrue)
        c(FPR = FPR(res$B[ix], Btrue[ix]), TPR = TPR(res$B[ix], Btrue[ix]) )
      }))
      auroc <- abs(AUROC(roc))
      name <- paste0("S_OU_", i, "_",
                     "p",p,"N",N, "d", d,
                     ".RData")
      save(file = paste0(path, "/", name), list = c("results", "roc", "exper", 
                                                    "auroc", "table", "p", 
                                                    "N", "d"))
    })  
  }
}

