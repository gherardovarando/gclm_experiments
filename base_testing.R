############ will remove it

library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

p <- 10
N <- 100000
d <- 2 / p

Btrue <- rStableMetzler(n = p, p = d, lower = FALSE, rfun = function(n) 
                         rnorm(n, sd = 1), 
                        rdiag = rnorm)

plot(graph_from_adjacency_matrix(sign(abs(t(Btrue)))))
Ctrue <- diag(p)
#Ctrue <- diag(runif(p, min = 0, max = 1))
#Ctrue <- diag(1 + runif(p, min = -0.5,2))
#Ctrue <- diag(p:1)

exper <- rOUinv(N, Btrue, C = Ctrue)

Sigmahat <- cor(exper$data)
#Sigmahat <- cov2cor(exper$Sigma)

########## path solutions
system.time(resllb <- llBpath(Sigmahat, eps = 1e-5, 
                               B0 = NULL,
                               lambdas = seq(0,2,length.out = 100),
                               maxIter = 5000, job = 11))

system.time(reslsb <- lsBpath(Sigmahat, eps = 1e-5, 
                  lambdas = seq(0,3,length.out = 100),
                  maxIter = 5000, job = 11))

reslasso <- lassoB(Sigmahat, C = diag(p))

resglasso <- glassoB(Sigmahat, 
                     lambda = seq(0, max(abs(Sigmahat)), length.out = 100))

evllb <- evaluatePathB(results = resllb, Btrue) 
evlsb <- evaluatePathB(results = reslsb, Btrue) 
evlasso <- evaluatePathB(results = reslasso, Btrue)
evglasso <- evaluatePathB(results = resglasso, Btrue)

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

deltaGraph( Btrue,resllb[[15]]$B, edge.arrow.size = 0.3, layout = layout_in_circle)
deltaGraph( Btrue,reslasso[[80]]$B, edge.arrow.size = 0.3, layout = layout_in_circle)

mllB(resllb[[which.max(evllb$confusion$f1)]]$B, S = (Sigmahat))
mll(solve(Sigmahat), Sigmahat)
################## Estimate C (B known)

p <- 10
Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
  rnorm(n, sd = 10), 
  rdiag = rnorm)
Ctrue <- diag(1 + runif(p, min = -0.5, max = 3))
exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))
B0 <- - 0.5 * diag(p) %*% solve(Sigmahat)

Sigmahat <- cov(exper$data)
Best <- proxgradllB(Sigmahat, B0, eps = 1e-14, lambda = 0.05)$B
system.time(out <- graddsllc(Sigmahat, B = Best, C = diag(p),C0 = diag(p),
                             maxIter = 10000,
                 eps = 1e-16, lambda = 0.005, alpha = 0.8, beta = 0.25))
plot(diag(Ctrue), diag(out$C))
abline(0,1)
diag(out$C)
max(abs(out$C - Ctrue))

#################### estimate B and C

p <- 20
N <- 10000
d <- 2 / p
Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
  rnorm(n, sd = 1), 
  rdiag = rnorm)
Ctrue <- diag(1 + runif(p, min = -0.5, max = 1))
#Ctrue <- diag(p:1)
#Ctrue <- diag(1:p/ sqrt(sum((1:p)^2)))
exper <- rOUinv(N, Btrue, C = Ctrue)

Sigmahat <- cov(exper$data)

C0 <- diag(p)
B0 <- - 0.5 * C0 %*% solve(Sigmahat)
deltaGraph(Btrue, Btrue, edge.arrow.size = 0.3, layout = layout_with_fr)
#B0 <- B0(p)
#B0 <- sign(abs(Btrue))
#diag(B0) <- -p
#B0 <- lowertriangB(Sigmahat)
out <- pnllbc(Sigma = Sigmahat, B = B0, C = C0, C0 = C0,
              eps = 1e-17, alpha = 0.5, beta = 0.25, maxIter = 1000, 
              intitr = 1,
              lambda = 0.15, lambdac = 0.0005, job = 11)
Best <- proxgradllB(Sigmahat, out$B, C = out$C, 
                    job = 10, eps = 1e-16, lambda = 0, maxIter = 1e3)$B
Cest <- graddsllc(Sigmahat, Best, out$C, eps = 1e-15, maxIter = 1000)$C
hamming(Btrue, out$B)
hamming(Btrue, Best)
deltaGraph(Btrue, out$B, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best, edge.arrow.size = 0.3, layout = layout_with_fr)
diag(out$C)
diag(Cest)
diag(Ctrue)
plot(diag(out$C), diag(Ctrue))
x <- diag(out$C / sqrt(sum(out$C^2)))
y <- diag(Ctrue/ sqrt(sum(Ctrue^2)))
plot(x,y)
abline(0,1,col= "red")

mllB(Best, Sigmahat, Cest)
mllB(out$B, Sigmahat, out$C)
mllB(Btrue, Sigmahat, Ctrue)

plot(out$B[Btrue != 0], Btrue[Btrue != 0])
abline(0,1,col= "red")

sum(Btrue!=0)
sum(Best!=0)
Best[Btrue == 0 & Best != 0]
