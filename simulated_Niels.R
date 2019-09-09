library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")



p <- 5
B <- matrix(nrow =  p, c(- 1, 1,  0,  0.2, 0,
                         -1, 0, 0.2, 0, 0,
                         0, 0, 0, -0.5, 0,
                         0, 0, 0, -1, 1,
                         0, 0, 1, 0, -1), byrow = T)


C <- diag(p)

exper <- rOUinv(1000, B)

B0 <- matrix(nrow = p, ncol = p, 1)
diag(B0) <- -p

Sigma <-  cov(exper$data)
B0 <- -0.5 * solve(Sigma)

res <- llBpath(Sigma, lambdas = seq(0,3, length.out = 100), job = 11,
               maxIter = 2000)

ev <- evaluatePathB(res, B)
plot(ev$roc)
AUROC(ev$roc)
ev$confusion

Best <- proxgradllB(Sigma = Sigma, B = Best, eps = 1e-15,
                   alpha = 0.5, maxIter = 1000,
                   lambda = 0,  job = 10)$B

hamming(B, Best)

deltaGraph(B, Best, edge.arrow.size = 0.3, layout = layout_with_fr)
mllB(Best, Sigma)
mllB(B, Sigma)
sum(B!=0)
sum(Best!=0)

