library(clggm)
library(igraph)
library(ggplot2)
library(glmnet)
source("functions/util.R")

p <- 10
N <- 10000
d <-  0.4
lower <- TRUE
M <- 50

lambda <- 0.1

res <- data.frame(t(replicate(M, {
  B <- rStableMetzler(n = p, p = d, lower = TRUE, rfun = function(n) 
    rnorm(n, sd = 10), 
    rdiag = rnorm)
  diag(B) <- diag(B) - 1
  C <- diag(p)
  gr <- graph_from_adjacency_matrix(abs(sign(t(B))),
                                    mode = "directed")
  ncomp <- components(gr)$no
  npar <- sum(B!=0)
  exper <- rOUinv(n = N, B = B, D = sqrt(C))
  Sigma <- cov(exper$data)
  #B0 <- lowertriangB(Sigma, C = C)
  #B0 <- -0.5 * C %*% solve(Sigma)
  B0 <- B0(p)
  t1 <- system.time(
      Best <-
        proxgradllB(Sigma = Sigma,
                    B = B0, eps = 1e-16, C = diag(p),
                    alpha = 0.5,
                    maxIter = 2000,
                    lambda = lambda,
                   job = 0)$B
    )[3]
  
  mllBest <- mllB(Best, Sigma)
  mllBt <- mll(solve(Sigma), Sigma)
  grest <- graph_from_adjacency_matrix(abs(sign(t(Best))))
  hamm <- hamming(B, Best)
  nparest <- sum(Best!=0)
  diffg <- union(difference(gr,grest), difference(grest,gr))
  deltaGraph(B, Best, edge.arrow.size = 0.3, layout = layout_with_fr,
             main = paste("lambda = ", lambda))
  return(c(
    ncomp = ncomp,
    npar = npar,
    nparest = nparest,
    hamm = hamm,
    mllB = mllBt,
    mllBest = mllBest,
    t1 = t1,
    densDiff = edge_density(diffg, loops = TRUE)
  ))
  })))

#dirname <- paste0("exp2", "-AlgLL","lambda", lambda,
#                  "p", p, "d", d, "N", N)

#if (lower){res <- data.frame(t(replicate(M, {
#  dirname <- paste0(dirname, "_lower")
#}

#dir.create(dirname, FALSE)

### minimum of the minus log likelihood versus estimated one
ggplot(data = res, aes(x = mllB, y = mllBest)) +
  geom_line(color = "lightblue") + geom_point() +
  ggtitle(paste("Log-likelihood comparison", "p =",p,", d =",d,", n =", N )) + xlab("- logLik full") +
  ylab("- logLik estimated")+
  geom_abline(slope = 1, intercept = 0, col = "red", size = 0.2)

#ggsave(paste0(dirname, "/loglik.pdf"))

###
ggplot(data = res) + geom_histogram(aes(x = hamm   ), bins = 10) +
  ggtitle("Distribution of Hamming distance")
#ggsave(paste0(dirname, "/hamming.pdf"))

ggplot(data = res) + geom_histogram(aes(x = hamm / (npar)   ), bins = 10) +
  ggtitle("Distribution of Hamming distance")
#ggsave(paste0(dirname, "/hammingfrac.pdf"))


ggplot(data = res) + geom_histogram(aes(x = hamm / (npar - p)   ), bins = 10) +
  ggtitle("Distribution of Hamming distance")
#ggsave(paste0(dirname, "/hammingfrachonest.pdf"))

ggplot(data = res, aes(x = npar, y = nparest)) +
  geom_line(color = "lightblue") + geom_point() +
  ggtitle(paste("Number of parameters", "p =",p,", d =",d,", n =", N )) + xlab("true") +
  ylab("estimated") +
  geom_abline(slope = 1, intercept = 0, col = "red", size = 0.2)

#ggsave(paste0(dirname, "/npars.pdf"))

#### time

ggplot(data = res) + geom_histogram(aes(x = t1.elapsed   ), bins = 10) +
  ggtitle("Distribution of elapsed time")
#ggsave(paste0(dirname, "/time.pdf"))

#save(res, file = paste0(dirname, "/res.RData"))

