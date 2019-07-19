library(clggm)
library(igraph)
source("functions/util.R")

p <- 10
N <- 1000
d <- 3 / p

Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
                         rnorm(n, sd = 10), 
                        rdiag = rnorm)
#diag(Btrue) <- diag(Btrue) - 1 

plot(graph_from_adjacency_matrix(sign(abs(t(Btrue)))))
Ctrue <- diag(p)
#Ctrue <- diag(1 + runif(p, min = -0.5,2))
#Ctrue <- diag(p:1)

exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))

Sigmahat <- cov(exper$data)


#######################################################################
lambda <- 0.05
C <- diag(p)
B0 <- - 0.5 * C %*% solve(Sigmahat)
#B0[abs(Sigmahat) < lambda^2] <- 0

### proximal gradient for -ll(B) + ||B||_1,off
### Best1 and Best2 should be the same (but sometimes no)
### Best1 is with the R implementation (gradient is computed in FORTRAN)
### Best2 is with pure fortran 
###  job = 0  ==  r = FALSE, h = FALSE
###  job = 10 == r = FALSE, h = TRUE
###  job = 11 == r = TRUE, h = TRUE
system.time(Best1 <- RproxgradllB(Sigma = Sigmahat, B = B0, C = C,
                                  eps = 1e-12,
                                  alpha = 1, beta = 0.5,
                                 maxIter = 1000, trace = 0, lambda = lambda,
                                  r = TRUE, h = FALSE))

system.time(Best2 <- proxgradllB(Sigmahat, B = B0, C = C,
                              eps = 1e-12,
                              alpha = 0.5,
                              maxIter = 1000, lambda = lambda,
                               job = 1)$B)

### proximal gradient for 0.5 * ||Sigma - S(B)||_2^2 + lambda * ||B||_1,off
system.time(Best3 <- proxgradlsB(Sigmahat, B = B0, C = C,
                                eps = 1e-5,
                                alpha = 0.5,
                                maxIter = 1000, lambda = lambda,
                                 job = 1)$B)

### coordinate descent for -ll(B) + ||B||_1,off
#system.time(Best4 <- proxcdllB(Sigmahat, B = B0, C = C,
#                                 eps = 1e-16,
#                                 alpha = 0.5,
#                                 maxIter = 1000, lambda = lambda,
#                                  job = 0)$B)
hamming(Best1, Best2)
hamming(Best1, Btrue)
hamming(Best2, Btrue)
hamming(Best3, Btrue)
#hamming(Best4, Btrue)
deltaGraph(Btrue, Best1, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best2, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best3, edge.arrow.size = 0.3, layout = layout_with_fr)
#deltaGraph(Btrue, Best4, edge.arrow.size = 0.3, layout = layout_with_fr)

########## path solutions
system.time(results <- llBpath(Sigmahat, eps = 1e-12,
                               maxIter = 1000, job = 11))
t(sapply(results, function(res) c(lambda = res$lambda, 
                                hamming = hamming(res$B, Btrue),
                                npar = sum(res$B!=0))))
