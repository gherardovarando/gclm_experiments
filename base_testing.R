library(clggm)
library(igraph)
source("functions/util.R")

p <- 20
N <- 1000000
d <- 2 / p

Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
                         rnorm(n, sd = 10), 
                        rdiag = rnorm)
diag(Btrue) <- diag(Btrue) - 1

plot(graph_from_adjacency_matrix(sign(abs(t(Btrue)))))
Ctrue <- diag(p)
#Ctrue <- diag(1 + runif(p, min = -0.5,2))
#Ctrue <- diag(p:1)

exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))

Sigmahat <- cov(exper$data)


#######################################################################
lambda <- 0.06
C <- diag(p)
B0 <- - 0.5 * C %*% solve(Sigmahat)
#B0[abs(Sigmahat) < lambda^2] <- 0

### proximal gradient for -ll(B) + ||B||_1,off
### Best1 and Best2 should be the same 
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
                               job = 11)$B)

### proximal gradient for 0.5 * ||Sigma - S(B)||_2^2 + lambda * ||B||_1,off
system.time(Best3 <- proxgradlsB(Sigmahat, B = B0, C = C,
                                eps = 1e-12,
                                alpha = 0.5,
                                maxIter = 1000, lambda = lambda,
                                 job = 0)$B)

### coordinate descent for -ll(B) + ||B||_1,off
system.time(Best4 <- proxcdllB(Sigmahat, B = B0, C = C,
                                 eps = 1e-16,
                                 alpha = 0.5,
                                 maxIter = 1000, lambda = lambda,
                                  job = 0)$B)
hamming(Best1, Best2)
hamming(Best1, Btrue)
hamming(Best2, Btrue)
hamming(Best3, Btrue)
hamming(Best4, Btrue)
deltaGraph(Btrue, Best1, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best2, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best3, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best4, edge.arrow.size = 0.3, layout = layout_with_fr)

########## path solutions
system.time(results <- llBpath(Sigmahat, eps = 1e-12,
                               maxIter = 1000, job = 11))

t(sapply(results, function(res) c(lambda = res$lambda, 
                                hamming = hamming(res$B, Btrue),
                                npar = sum(res$B!=0))))


################## Estimate C (B known)

p <- 10
Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
  rnorm(n, sd = 10), 
  rdiag = rnorm)
Ctrue <- diag(1 + runif(p, min = -0.5, max = 3))
exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))
B0 <- - 0.5 * Ctrue %*% solve(Sigmahat)

Sigmahat <- cov(exper$data)

system.time(out <- graddsllc(Sigmahat, B = Btrue, C = diag(p),C0 = diag(p),
                             maxIter = 10000,
                 eps = 1e-16, lambda = 0.0005, alpha = 0.5, beta = 0.2))
plot(diag(Ctrue), diag(out$C))
diag(out$C)
max(abs(out$C - Ctrue))

#################### estimate B and C

p <- 10
N <- 10000
d <- 2 / p
Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
  rnorm(n, sd = 10), 
  rdiag = rnorm)
Ctrue <- diag(1 + runif(p, min = -0.5, max = 1))
#Ctrue <- diag(1:p/ sqrt(sum((1:p)^2)))
exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))

Sigmahat <- cov(exper$data)

C0 <- diag(p)
B0 <- - 0.5 * C0 %*% solve(Sigmahat)
#B0 <- B0(p)
out <- pnllbc(Sigma = Sigmahat, B = B0, C = C0, C0 = C0,
              eps = 1e-8, alpha = 0.5, beta = 0.25, maxIter = 100, 
              intitr = 1,
              lambda = 0.35, lambdac = 0.005, job = 0)
hamming(Btrue, out$B)
deltaGraph(Btrue, out$B, edge.arrow.size = 0.3, layout = layout_with_fr)
out$C
Ctrue
plot(diag(out$C), diag(Ctrue))
diag(out$C - Ctrue)
x <- diag(out$C / sqrt(sum(out$C^2)))
y <- diag(Ctrue/ sqrt(sum(Ctrue^2)))
plot(x,y)
abline(0,1,col= "red")

