library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

N <- 100000
p <- 10
B <- matrix(nrow = p, ncol = p, 0 )
diag(B) <-  - ((p+1):(2)) / 5

for (i in 2:p){
  B[i-1, i] <- 1
}
B[i, 1] <- 1

C <- diag(runif(p))
C <- diag(p)
exper <- rOUinv(N, B, D = sqrt(C))

Sigmahat <- cor(exper$data)

results <- llBpath(Sigmahat, eps = 1e-15, C = C,
                   maxIter = 5000, job = 11, 
                   lambdas = seq(0,3,length.out = 100))
roc <- t(sapply(results, function(res){
  ix <- lower.tri(B) | upper.tri(B)
  c(FPR = FPR(res$B[ix], B[ix]), TPR = TPR(res$B[ix], B[ix]) )
}))
auroc <- abs(AUROC(roc))
plot(roc)

ggplot(data=as.data.frame(roc), aes(x= FPR, y= TPR)) +
  geom_path()+
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") + 
  coord_fixed() + 
  ggtitle(paste("ROC"), subtitle = paste0("p=", p, ", N=", N))+
  ggsave(paste0("plot/", "cycle_N", N, "p", p, "roc_envelope.pdf"))

t(sapply(results, function(res){
  return(c(lambda = res$lambda, hamming = hamming(B, res$B)))
}))

Sigmahat <- exper$Sigma
B0 <- -0.5 * solve(Sigmahat)
Best <- proxgradllB(Sigmahat, B0, C =  C, eps = 1e-12, maxIter = 5000, 
                    lambda = 0.5, job = 1)$B
deltaGraph(B, Best, layout = layout_in_circle, edge.arrow.size = 0.3)


res <- pnllbc(Sigmahat, B0, eps = 1e-15, maxIter = 1000, intitr = 10, 
       lambda = 0.025, lambdac = 0.0005, job = 1)

deltaGraph(B, res$B, layout = layout_in_circle, edge.arrow.size = 0.3)
diag(res$C)
diag(C)
plot(diag(C), diag(res$C))
abline(0,1)
