library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

N <- 10000
p <- 20
B <- matrix(nrow = p, ncol = p, 0 )
diag(B) <-  - 2*(1:p)

for (i in 2:p){
  B[i-1, i] <- 1
}
B[i, 1] <- 1
B[p/2, 1] <- 1
C <- diag(runif(p))
exper <- rOUinv(N, B, C = C)

Sigmahat <- cor(exper$data[1:5000,])

results <- llBpath(Sigmahat, lambdas = seq(0,3,length.out = 100))
