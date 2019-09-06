library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

p <- 10
N <- 1000
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

system.time(Best <- proxgradllB(Sigmahat, B = B0, C = C,
                              eps = 1e-12,
                              alpha = 0.5,
                              maxIter = 1000, lambda = lambda,
                               job = 0)$B)

### proximal gradient for 0.5 * ||Sigma - S(B)||_2^2 + lambda * ||B||_1,off
system.time(Best2 <- proxgradlsB(Sigmahat, B = B0, C = C,
                                eps = 1e-12,
                                alpha = 0.5,
                                maxIter = 1000, lambda = lambda,
                                 job = 10)$B)

### coordinate descent for -ll(B) + ||B||_1,off
system.time(Best3 <- proxcdllB(Sigmahat, B = B0, C = C,
                                 eps = 1e-16,
                                 alpha = 0.5,
                                 maxIter = 1000, lambda = lambda,
                                  job = 0)$B)

hamming(Best, Btrue)
hamming(Best2, Btrue)
hamming(Best3, Btrue)
deltaGraph(Btrue, Best, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best2, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best3, edge.arrow.size = 0.3, layout = layout_with_fr)

mllB(Best, Sigmahat, C = Ctrue)
mllB(Btrue, Sigmahat, C = Ctrue)

### refine with maximum likelihood
BestI <- proxgradllB(Sigmahat, Best, eps = 1e-16, job = 10)$B
### correct structure
Bs <- sign(abs(Btrue))
diag(Bs) <- -p
BestII <- proxgradllB(Sigmahat, Bs, job = 10, eps = 1e-16)$B
BestI
BestII
D <- diag(runif(p))
BestIII <- proxgradllB(D %*% Sigmahat %*% D, Bs,C = D %*% D, 
                       job = 10, eps = 1e-16)$B
B4 <- proxgradllB(D %*% Sigmahat %*% D, Bs,C = diag(p), 
                  job = 10, eps = 1e-16, maxIter = 20000)$B

xx <- diag(Btrue) / diag(B4)
plot(xx, diag(D)^2)
abline(0,1)
########## path solutions
system.time(results <- llBpath(Sigmahat, eps = 1e-12,
                               maxIter = 1000, job = 10))

t(sapply(results, function(res) c(lambda = res$lambda, 
                                hamming = hamming(res$B, Btrue),
                                npar = sum(res$B!=0))))

resultspath <- llBpath(cov2cor(Sigmahat), eps = 1e-12, 
                       maxIter = 5000, job = 10, 
                       lambdas = seq(0,2,length.out = 100))

t(sapply(resultspath, function(res) c(lambda = res$lambda, 
                                  hamming = hamming(res$B, Btrue),
                                  npar = sum(res$B!=0))))

roc <- t(sapply(resultspath, function(res){
  ix <- lower.tri(Btrue) | upper.tri(Btrue)
  c(FPR = FPR(res$B[ix], Btrue[ix]), TPR = TPR(res$B[ix], Btrue[ix]) )
}))
plot(roc, type = "l")
AUROC(roc)
ggplot(data=data.frame(roc), aes(x= FPR, y= TPR)) +
  geom_path()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") 
  #+
  #coord_fixed()
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

p <- 10
N <- 10000
d <- 2 / p
Btrue <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
  rnorm(n, sd = 10), 
  rdiag = rnorm)
Ctrue <- diag(1 + runif(p, min = -0.5, max = 1))
#Ctrue <- diag(p:1)
#Ctrue <- diag(1:p/ sqrt(sum((1:p)^2)))
exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))

Sigmahat <- cov(exper$data)

C0 <- diag(p)
B0 <- - 0.5 * C0 %*% solve(Sigmahat)
deltaGraph(Btrue, Btrue)
#B0 <- B0(p)
#B0 <- sign(abs(Btrue))
#diag(B0) <- -p
#B0 <- lowertriangB(Sigmahat)
out <- pnllbc(Sigma = Sigmahat, B = B0, C = C0, C0 = C0,
              eps = 1e-12, alpha = 0.5, beta = 0.25, maxIter = 1000, 
              intitr = 1000,
              lambda = 0.2, lambdac = 0, job = 10)
Best <- proxgradllB(Sigmahat, out$B, C = out$C, 
                    job = 10, eps = 1e-16, lambda = 0, maxIter = 1e5)$B
Cest <- graddsllc(Sigmahat, Best, out$C, eps = 1e-15, maxIter = 1000)$C
hamming(Btrue, out$B)
hamming(Btrue, Best)
deltaGraph(Btrue, out$B, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best, edge.arrow.size = 0.3, layout = layout_with_fr)
diag(out$C)
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


