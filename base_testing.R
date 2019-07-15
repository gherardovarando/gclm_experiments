library(clggm)
library(igraph)
source("functions/util.R")

p <- 50
N <- 1000
d <- 0.07

Btrue <- rStableMetzler(p, d = d, lower = TRUE, rfun = runif)

plot(graph_from_adjacency_matrix(sign(abs(t(Btrue)))))
Ctrue <- diag(p)
#Ctrue <- diag(1 + runif(p, min = -0.5,2))
#Ctrue <- diag(p:1)

exper <- rOUinv(N, Btrue, D = sqrt(Ctrue))

Sigmahat <- cov(exper$data)


#######################################################################
lambda <- 0.05
Cest <- diag(p)
Best <- - 0.5 * Cest %*% solve(Sigmahat)
Best[abs(Sigmahat) < lambda^2] <- 0
#Cest <- estimateCLL(Sigmahat, B = Best, C = Cest,
#                                  eps = 1e-10, alpha = 1, beta = 0.5,
#                                  maxIter = 100, trace = 2, lambda = 0.1,
#                                  t0 = 1)

system.time(Best2 <- proxgradB(Sigma = Sigmahat, B = Best, C = Cest,
                                  eps = 1e-8,
                                  alpha = 1, beta = 0.5,
                                 maxIter = 1000, trace = 0, lambda = lambda,
                                  r = FALSE, h = FALSE))


system.time(Best <- fproxgradB(Sigmahat, B = Best, C = Cest,
                              eps = 1e-8,
                              alpha = 0.5,
                              maxIter = 3, lambda = lambda,
                              all = FALSE, job = 11))

mllB(Btrue, Sigmahat, C = exper$C) + lambda * (sum(abs(Btrue)) - sum(abs(diag(Btrue))))
mllB(Best2, S =  Sigmahat, C = Cest) + lambda * (sum(abs(Best2)) - sum(abs(diag(Best2))))
mllB(Best, S =  Sigmahat, C = Cest) + lambda * (sum(abs(Best)) - sum(abs(diag(Best))))

hamming(Best2, Best)
hamming(Best2, Btrue)
hamming(Btrue, Best)
deltaGraph(Btrue, Best2, edge.arrow.size = 0.3, layout = layout_with_fr)
deltaGraph(Btrue, Best, edge.arrow.size = 0.3, layout = layout_with_fr)


###############################################################
lambda <- 0.03
Cest <- diag(p)
Best <- -0.5 * Cest %*% solve(Sigmahat)

#Best <- lowertriangB(Sigmahat, C= exper$C)
for (i in 1:100){
  S <- cov(exper$data[sample(1:N, size = N / 2, replace = FALSE),])

  Best <-  proxgradB(S, B = Best, C = Cest,
                                   eps = 1e-8,
                                   alpha = 1, beta = 0.5,
                                   maxIter = 10, trace = 1, lambda = lambda,
                                   r = TRUE, h = FALSE)
  #Cest <- estimateCLL(S, B = Best, C = Cest,
  #                                 eps = 1e-10, alpha = 1, beta = 0.5,
  #                                 maxIter = 100, trace = 1, lambda = 0.1,
  #                                 t0 = 1)

}
mllB(Btrue, C = exper$C, S = Sigmahat) + lambda * (sum(abs(Btrue)) - sum(abs(diag(Btrue))))
mllB(Best, S =  Sigmahat, C = Cest)+ lambda * (sum(abs(Best)) - sum(abs(diag(Best))))

hamming(Btrue, Best)
deltaGraph(Btrue, Best, edge.arrow.size = 0.3, layout = layout_with_fr)


Best <- sign(abs(Best))
Best[Best == 0 ] <- 1e-6
Best[Best > 1] <- 1
AUROC(Best, Btrue, 100)
roc <- ROC(Best, Btrue)
ggplot(data=data.frame(roc), aes(x= FPR, y= TPR)) +
  geom_line()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  coord_fixed()


#####################################################################
### bootstrapping for confidence
lambda <- 0.01
Bp <- matrix(nrow = p, ncol = p, 0)

for (i in 1:100){
  S <- cov(exper$data[sample(1:N, size = N / 2, replace = FALSE),])
  Best <- -0.5 * exper$C %*% solve(S)
  Bp <- Bp +  sign(abs(fproxgradB(Sigma = S, B = Best, C = exper$C,
                                  lambda = lambda,
                      eps = 1e-8,alpha =0.5, maxIter = 100, job = 0, 
                      all = FALSE)))
}

Bp <- Bp / 100

deltaGraph(Btrue, Bp, edge.arrow.size = 0.3, layout = layout_with_fr)

Bp[Bp == 0] <- 1e-6

ix <- lower.tri(Bp) | upper.tri(Bp)
AUROC(Bp[ix], Btrue[ix])

roc <- ROC(Bp[ix], Btrue[ix])
ggplot(data=data.frame(roc), aes(x= FPR, y= TPR)) +
  geom_line()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  coord_fixed()


############################################## Cest
x <- diag(Cest) / sqrt(sum(diag(Cest^2)))
y <- diag(exper$C) / sqrt(sum(diag(exper$C^2)))

sqrt(sum( (x - y)^2 ))
acos(sum(x * y))
plot(x,y)
abline(0,1, col = "red")
x
y

