library(clggm)
library(ggplot2)
library(igraph)
library(glmnet)
source("functions/util.R")

##################### read data and true graph 
D <- read.csv("data/dreamhpn/insilico/CSV/insilico.csv")
trueGraph <-
  as.matrix(read.csv("data/dreamhpn/SC1B_GS/SC1B_GS/trueGraph.csv",header = FALSE,as.is = TRUE   ))
p <- 20
diag(trueGraph) <- 1
colnames(trueGraph) <- rownames(trueGraph) <- NULL


###################### estimate gradph from all data
S <- cor(D[, -(1:4)])
B0 <- - 0.5 * diag(p) %*% solve(S)
C0 <- diag(p)

Best <- proxgradllB(Sigma = S,B = B0, C = C0, eps = 1e-10, 
                   alpha = 0.5, maxIter = 1000, lambda = 0.5, job = 0)$B

deltaGraph(t(trueGraph), Best,
           edge.arrow.size = 0.3, layout = layout_in_circle)

hamming(t(trueGraph), Best)

results <- llBpath(S, job = 11)
B <- t(trueGraph)
t(sapply(results, function(res) c(lambda = res$lambda, 
                                  npar = sum(res$B!=0),
                                  fp = sum(res$B!=0 & B==0),
                                  tp = sum(res$B!=0 & B!=0) ,
                                  fn = sum(res$B==0 & B!=0),
                                  tn = sum(res$B==0 & B==0),
                                  errs = sum(res$B!=0 & B==0) + 
                                    sum(res$B==0 & B!=0))))
################# 
S <- cor(D[, -(1:4)])
Best <- - 0.5 * diag(p) %*% solve(S)
C0 <- diag(p)
lambda <- 0.05
N <- nrow(D)
for (i in 1:100){
  S <- cor(D[sample(1:N, size = N, replace = TRUE), - (1:4)])
  Best <-  proxgradllB(Sigma = S,B = Best, C = C0, eps = 1e-10, 
                      alpha = 0.5, maxIter = 100, lambda = lambda,
                      job = 1)$B
  print(i)
}

deltaGraph(t(trueGraph), Best,
           edge.arrow.size = 0.3, layout = layout_in_circle)
hamming(Best, t(trueGraph))

Bp <- sign(abs(Best))
Bp[Bp == 0] <- 1e-9
Bp[Bp > 1] <- 1
AUROC(Bp , t(trueGraph))
roc1 <- ROC(Bp , t(trueGraph))

ggplot(data=data.frame(roc1), aes(x= FPR, y= TPR)) +
  geom_line()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  coord_fixed()

##################################

# library(parallel)
# cl <- makeCluster(1)
# clusterExport(cl, "D")
# clusterExport(cl, "p")
# clusterCall(cl, function() library(clggm) )
res <- sapply(X = 1:100, FUN = function(x){
  S <- cor(D[sample(1:nrow(D), replace = TRUE, size = nrow(D)  ), 
             -(1:4)])
  B0 <- - 0.5 * diag(p) %*% solve(S)
  C <- diag(p)
  Best <- proxgradllB(Sigma = S, C = C, B = B0,
                        lambda = 0.03, 
                        alpha = 0.5,eps = 1e-8, 
                        maxIter = 1000, job = 1)$B
  print(x)
  return(Best)
})

AA <- apply(res, MARGIN = 1, FUN =  function(x) mean((abs(x))))

AA <- matrix(nrow = 20, ncol = 20, AA)

deltaGraph(t(trueGraph), AA,
           edge.arrow.size = 0.3, layout = layout_in_circle)

Bp <- AA
Bp[Bp == 0] <- 1e-9
Bp[Bp > 1] <- 1
ix <- lower.tri(Bp) | upper.tri(Bp)
AUROC(Bp[ix] , t(trueGraph)[ix], 300)
roc1 <- ROC(Bp[ix] , t(trueGraph)[ix])

ggplot(data=data.frame(roc1), aes(x= FPR, y= TPR)) +
  geom_line()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  coord_fixed()
################################################################

p <- 20
C <- diag(p)
lambda <- 0.15
Btot <- matrix(nrow = p, ncol = p, 0)
for (stimulus in levels(D$Stimulus)){
  for (inh in levels(D$Inhibitor)){
    if (nrow(D[D$Stimulus == stimulus & D$Inhibitor == inh, -(1:4)]) > 1){
      S <- cor(D[D$Stimulus == stimulus & D$Inhibitor == inh, -(1:4)])
      Best <- proxgradllB(Sigma = S, C = C, B = B0,
                          lambda =  lambda, 
                          alpha = 0.5,eps = 1e-8, 
                          maxIter = 1000, job = 1)$B
      deltaGraph(t(trueGraph), Best,
                 edge.arrow.size = 0.3, layout = layout_with_fr)
      #Btot <- Btot + sign(abs(Best))
      Btot <- Btot + abs(Best) / sum(diag(clyap(Best, c(diag(p)))))
    }
    
  }
}


Btot2 <- (Btot)
Btot2[Btot2 <= 0] <- 1e-6
Btot2[Btot2 > 1] <- 1
ix <- lower.tri(Btot) | upper.tri(Btot)
AUROC(Btot2[ix], t(trueGraph)[ix], 100)
roc2 <- ROC(Btot2[ix], t(trueGraph)[ix])

ggplot(data=data.frame(roc2), aes(x= FPR, y= TPR)) +
  geom_line()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  coord_fixed()

#######################################################################
p <- 20
Btot <- matrix(nrow = p, ncol = p, 0)
  for (stimulus in levels(D$Stimulus)){
  for (inh in levels(D$Inhibitor)){
    if (nrow(D[D$Stimulus == stimulus & D$Inhibitor == inh, -(1:4)]) > 1){
      S <- cor(D[D$Stimulus == stimulus & D$Inhibitor == inh, -(1:4)])
      MM <- matrix(nrow = p, ncol = p, 1:(p^2))
      TT  <- diag(p ^ 2)[c(t(MM)),]
      AA <- S %x% diag(p) + ( (diag(p) %x% S) %*% TT)
      Best <- matrix(nrow = p, ncol = p,
                   glmnet(AA, y = - c(diag(p)),
                          lambda = 0.0008,
                          intercept = FALSE,
                          standardize = TRUE)$beta[,1])
      deltaGraph(t(trueGraph), Best,
                 edge.arrow.size = 0.3, layout = layout_with_fr)
      #Btot <- Btot + sign(abs(Best))
      Btot <- Btot + abs(Best) / sum(diag(clyap(Best, c(diag(p)))))
    }

  }
}


Btot2 <- (Btot)
Btot2[Btot2 <= 0] <- 1e-6
Btot2[Btot2 > 1] <- 1
ix <- lower.tri(Btot) | upper.tri(Btot)
AUROC(Btot2[ix], t(trueGraph)[ix], 100)
roc2 <- ROC(Btot2[ix], t(trueGraph)[ix])

ggplot(data=data.frame(roc2), aes(x= FPR, y= TPR)) +
  geom_line()+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  coord_fixed()

