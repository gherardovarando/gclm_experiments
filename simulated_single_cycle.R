library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

N <- 100000
p <- 10
B <- matrix(nrow = p, ncol = p, 0 )
diag(B) <-  - 3

for (i in 2:p){
  B[i-1, i] <- 1
}
B[i, 1] <- 1

C <- diag(p)
exper <- rOUinv(N, B, C = C)

Sigmahat <- cor(exper$data)

resllb <- llBpath(Sigmahat, eps = 1e-6, B0 = B0(p),
                   maxIter = 5000, job = 11, 
                   lambdas = seq(0,3,length.out = 100))

reslsb <- lsBpath(Sigmahat, eps = 1e-6, 
                              lambdas = seq(0,0.5,length.out = 100),
                              maxIter = 5000, job = 11)

reslasso <- lassoB(Sigmahat, C = diag(p))

resglasso <- glassoB(Sigmahat)

evllb <- evaluatePathB(results = resllb, B) 
evlsb <- evaluatePathB(results = reslsb, B) 
evlasso <- evaluatePathB(results = reslasso, B)
evglasso <- evaluatePathB(results = resglasso, B)

AUROC(evllb$roc)
AUROC(evlsb$roc)
AUROC(evlasso$roc)
AUROC(evglasso$roc)

plotROC(evllb$roc)
plotROC(evlsb$roc)
plotROC(evlasso$roc)

plotROCS(list(logLik_l1 = evllb, 
              ls_l1 = evlsb, 
              LASSO = evlasso,
              GLASSO = evglasso))

plotPR(evllb$confusion)
plotPR(evlsb$confusion)
plotPR(evlasso$confusion)
plotPR(evglasso$confusion)

max(evllb$confusion$f1)
max(evlsb$confusion$f1)
max(evlasso$confusion$f1)
max(evglasso$confusion$f1)

max(evllb$confusion$bacc)
max(evlsb$confusion$bacc)
max(evlasso$confusion$bacc)
max(evglasso$confusion$bacc)

min(evllb$confusion$errs)
min(evlsb$confusion$errs)
min(evlasso$confusion$errs)
min(evglasso$confusion$errs)

Best <- proxgradllB(Sigma = Sigmahat, 
                    B = B, lambda = 0,
                    job =10, eps = 1e-16, maxIter = 15000)$B
deltaGraph( B, Best, edge.arrow.size = 0.2, layout = layout_with_fr)
mllB(Best, S = Sigmahat)
mll(solve(Sigmahat), Sigmahat)

