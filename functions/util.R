#' Difference graph
#'
#' @param B1 matrix
#' @param B2 matrix
#' @param legend logical 
#' @param ... parameters to be passed to the igraph plotting
#'
#' @return a graph with colored edges where
#'         in red are plotted the edges in \code{B1} and not in
#'         \code{B2} and in blue the edges in \code{B2} and not in
#'         \code{B1}.
#' @export
deltaGraph <- function(B1, B2, legend = FALSE,  ...){
  gr1 <- igraph::graph_from_adjacency_matrix(abs(sign(t(B1))))
  gr2 <- igraph::graph_from_adjacency_matrix(abs(sign(t(B2))))
  d1 <-  igraph::set_edge_attr(gr1 - gr2, name = "color", index = , value = "red" )
  d1 <- igraph::add_edges(d1, edges = c(t(igraph::get.edgelist(gr2 - gr1))),
                  attr = list(color = "blue"))
  ed <- c(t(igraph::get.edgelist(igraph::intersection(gr1, gr2))))
  d1 <- igraph::add_edges(d1, edges = ed,
                  attr = list(color = "black"))
  igraph::plot.igraph(d1, ...)
  if (legend){
    legend("left", legend = c("correct", "missing", "wrong"),
           col = c("black", "red", "blue"), lty = 1, lwd = 2 ,cex = 0.5)
  }
  return(invisible(d1))
}


#' Hamming distance between adjacency matrices
#'
#' @param A matrix
#' @param B matrix
#'
#' @return The hamming distance between the directed graphs induced by
#'  \code{A} and \code{B}.
#' @export
hamming <- function(A, B){
  sum( (A != 0) != (B != 0) )
}


#' Minus Log-Likelihood for Gaussian model
#'
#' @param P the inverse of the covariance matrix
#' @param S the empirical covariance
#'
#' @importFrom stats var
#' @export
mll <- function(P, S){
 -determinant(P, logarithm = TRUE)$modulus + sum(S * P)
}

#' Minus log-likelihood for the invariant distribution of OU
#'
#' @param B coefficent matrix
#' @param S covariance matrix
#' @param C noise matrix
#'
#' @export
mllB <- function(B, S, C = diag(nrow(B))){
  P <- solve(clyap(B, C))
  mll(P, S)
}


#' Generate data from invariant distribution of O-U process
#'
#' Sample observations from the invariant distribution of an
#' Ornstein-Uhlenbeck process:  \deqn{dX_t = B(X_t - \mu)dt + DdW_t}
#'
#' @param n number of sample
#' @param B square matrix, coefficients of the O-U equation
#' @param D noise matrix
#' @param mean invariant mean
#'
#' @return  A list with components:
#'   -  \code{sample} the produced sample
#'   -  \code{Sigma} the covariance matrix of the invariant distribution
#'   -  \code{C} the  C matrix obtained as DD'
#'   -  \code{B} the coefficient square matrix
#'   -  \code{mean} the mean vector of the invariant distribution
#'
#'
#' @importFrom MASS mvrnorm
#'
#' @export
rOUinv <- function(n = 1, B,
                C = diag(nrow = nrow(B), ncol = ncol(B)),
                mean = rep(0, nrow(B)) ){
  Sigma <- clyap(B, C)
  S <- MASS::mvrnorm(n = n, Sigma = Sigma, mu = mean)
  return(list(data = S, Sigma = Sigma, C = C, B = B, mean = mean))
}




#' Generate a random skew-symmetric matrix
#'
#' @param p integer the dimension of the desired matrix
#' @param rfun a function with signature \code{function(n)} that
#' generates \code{n} numbers (e.g. \code{rnorm})
#'
#' @return a skew-symmetric matrix
#'
#' @importFrom stats rnorm
#' @export
rSkewSymm <- function(p, rfun = rnorm){
  M <- matrix(nrow = p, ncol = p, data = rfun(p * p))
  return(M - t(M) )
}

#' Generate stable Metzler matrix
#'
#' @param n integer the dimension of the matrix
#' @param p probability of non zero entries
#' @param lower logical, should the matrix be lower triangular
#' @param rfun the random number generator for the matrix entries
#' @param rdiag the random number generator for the diagonal added error
#'
#' @return a square matrix with positive off-diagonal entries (Metzler)
#' and eigenvalues with strict negative real part (stable)
#'
#' @importFrom stats rnorm
#' @export
rStableMetzler <- function(n = 1, p = 1, lower = FALSE,
                           rfun = rnorm, rdiag = rnorm){
  D <- abs(rfun(n * n)) * sample(c(0,1), replace = TRUE, size = n * n,
                             prob = c(1 - p, p))
  A <- matrix(nrow = n, ncol = n, data = D)
  diag(A) <- 0
  if (lower){
    A[upper.tri(A)] <- 0
  }
  diag(A) <- -rowSums(A) - abs(rdiag(n))
  return(A)
}


#' Compute the ROC curve 
#' 
#' @param x an array with confidence levels between 0 and 1
#' @param y the true target values (0 and >0)
#' @param n number of points 
#' @export
ROC <- function(x, y, n = 100){
 pp <- seq(from =1, to = 0, length.out = n)
 D <- sapply(pp, FUN = function(p){
   x[x <= p] <- 0
   c(FPR = FPR(x, y), TPR = TPR(x, y))
 })
 return(t(D))
}

#'@rdname ROC
#'@export
AUROC <- function(x, y = NULL, n = 100){
 if (!is.null(y)){
   D <- ROC(x, y, n)   
 }else{
   D <- x
   n <- nrow(D)
 }
 temp <- 0
 for (i in 1:(n - 1)){
   ## trapezoidal quadrature rule
   temp <- temp + (D[i + 1,2] + D[i, 2]) * ( D[i + 1,1] - D[i,1]) / 2
 }
 return(temp)
}



#'@rdname ROC
#'@export
FPR <- function(x, y){
  if (sum(y == 0) == 0){
    return(1)
  }
  (sum(x != 0 & y == 0)) / (sum(y == 0))
}


#'@rdname ROC
#'@export
TPR <- function(x, y){
  if (sum(y != 0) == 0){
    return(1)
  }
   sum(x != 0 & y != 0) / sum(y !=0) 
}


#' Generate a naive stable matrix 
#' 
#' @param p dimension of the matrix
#' @return a stable matrix with off-diagonal entries equal to 1 and 
#' diagonal entries equal to \code{-p}
#' @export
B0 <- function(p){
  M <- matrix(nrow = p, ncol = p, 1)
  diag(M) <- -p
  return(M)
}


evaluatePathB <- function(results, B){
  ix <- lower.tri(B) | upper.tri(B)
  conf <- as.data.frame(t(sapply(results, function(res) c(lambda = res$lambda, 
                                    npar = sum(res$B[ix]!=0),
                                    fp = sum(res$B[ix] !=0 & B[ix] ==0),
                                    tp = sum(res$B[ix] !=0 & B[ix] !=0) ,
                                    fn = sum(res$B[ix] ==0 & B[ix] !=0),
                                    tn = sum(res$B[ix] ==0 & B[ix] ==0),
                                    errs = sum(res$B[ix] !=0 & B[ix] ==0) + 
                                      sum(res$B[ix] == 0 & B[ix] !=0)))))
  conf$tpr <- conf$tp / (conf$tp + conf$fn)
  conf$tnr <- conf$tn / (conf$tn + conf$fp)
  conf$bacc <- (conf$tpr + conf$tnr) / 2
  conf$precision <- conf$tp / (conf$tp + conf$fp)
  conf$precision[is.nan(conf$precision)] <- 1
  conf$recall <- conf$tp / (conf$tp + conf$fn)
  conf$recall[is.nan(conf$recall)] <- 1
  conf$f1 <- 2 * (conf$precision * conf$recall) / (conf$precision + conf$recall)
  conf$f1[is.nan(conf$f1)] <- 0
  roc <- t(sapply(results, function(res){
    c(FPR = FPR(res$B[ix], B[ix]), TPR = TPR(res$B[ix], B[ix]) )
  }))
  roc <- roc[dim(roc)[1]:1,]
  roc <- rbind(c(0,0), roc, c(1,1))
  return(list(roc = roc, confusion = conf))
}


library(lars)
library(glmnet)
lassoB <- function(Sigma, C = diag(nrow(Sigma)), lambda = NULL){
  p <- nrow(Sigma)
  MM <- matrix(nrow = p, ncol = p, 1:(p^2))
  TT  <- diag(p ^ 2)[c(t(MM)),]
  AA <- Sigma %x% diag(p) + ( (diag(p) %x% Sigma) %*% TT)
#  tmp <- lars(AA, y = - c(C),
#               type = "lasso",
#                intercept = FALSE,
#                normalize = FALSE )
  if (is.null(lambda)){
    tmp <- glmnet(AA, y = -c(C),
                  intercept = FALSE,
                  standardize = FALSE,
                  nlambda = 3,
                  penalty.factor = 1 - diag(p))
    lambdamax <- max(tmp$lambda)
    lambda <- seq(lambdamax, 0, length.out =  100)
  }

  tmp <- glmnet(AA, y = -c(C),
                intercept = FALSE,
                standardize = FALSE,
                nlambda = 100,
                lambda = lambda,
                lambda.min.ratio = 1e-8,
                penalty.factor = 1 - diag(p),
                maxit = 1e+8)
### ATTENTION if use glmnet change to tmp$beta[,i]:
### if use lars change to tmp$beta[i,]
    obj <- lapply(length(tmp$lambda):1, function(i){
    list(B = matrix(nrow =p, ncol = p, data = tmp$beta[,i]), 
         lambda = tmp$lambda[i])
  })
    attr(obj, "jerr") <- tmp$jerr
    attr(obj, "nulldev") <- tmp$nulldev
    return(obj)
}


library(glasso)
glassoB <- function(Sigma, lambda = NULL){
  gpath <- glassopath(Sigma, rholist = lambda, trace = FALSE)
  return(lapply(1:length(gpath$rholist), FUN = function(i){
    list(B = gpath$wi[,,i], lambda = gpath$rholist[i])
  }))
}


library(ggplot2)

plotROC <- function(roc){
  ggplot(data=as.data.frame(roc), aes(x= FPR, y= TPR)) +
    geom_path()+
    geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") + 
    coord_fixed()
}


plotPR <- function(x){
  ggplot(data=x, aes(y= precision, x= recall)) +
    geom_path()+
    geom_abline(intercept = 1, slope = -1, col = "gray", linetype = "dashed") + 
    coord_fixed()
}

plotROCS <- function(x, title = NULL){
  nms <- names(x)
  df <- as.data.frame(x[[1]]$roc)
  df$algorithm <- nms[1]
  for (j in 2:length(x)){
    dtmp <- as.data.frame(x[[j]]$roc)
    dtmp$algorithm <- nms[j]
    df <- rbind(df, dtmp)
  }
  ggplot(data=df, aes(x= FPR, y= TPR, group = algorithm, color = algorithm)) +
    geom_path()+
    geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") + 
    coord_fixed() + ggtitle(title) 
}


igraph.to.tikz <- function (graph, layout) {
  ## Here we get the matrix layout
  if (class(layout) == "function")
    layout <- layout(graph)
  
  layout <- layout / max(abs(layout))
  
  ##TikZ initialisation and default options (see pgf/TikZ manual)
  cat("\\tikzset{\n")
  cat("\tnode/.style={circle,inner sep=1mm,minimum size=0.8cm,draw,very thick,black,fill=red!20,text=black},\n")
  cat("\tnondirectional/.style={very thick,black},\n")
  cat("\tunidirectional/.style={nondirectional,shorten >=2pt,-stealth},\n")
  cat("\tbidirectional/.style={unidirectional,bend right=10}\n")
  cat("}\n")
  cat("\n")
  
  ##Size
  cat("\\begin{tikzpicture}[scale=5]\n")
  
  for (i in 1:length(V(graph))) {
    vertex <- V(graph)[i]
    label <- V(graph)[vertex]$name
    if (is.null(label))
      label <- ""
    
    ##drawing vertices
    cat (sprintf ("\t\\node [node] (v%d) at (%f, %f)\t{%s};\n", vertex, layout[i,1], layout[i,2], label))
  }
  cat("\n")
  
  adj = get.adjacency(graph)
  bidirectional = adj & t(adj)
  
  if (!is.directed(graph)) ##undirected case
    for (line in 1:nrow(adj)) {
      for (col in line:ncol(adj)) {
        if (adj[line,col]&col>line) {
          cat (sprintf ("\t\\path [nondirectional] (v%d) edge (v%d);\n", line, col)) ##edges drawing
        }
      }
    }
  else ##directed case
    for (line in 1:nrow(adj)) {
      for (col in 1:ncol(adj)) {
        if (bidirectional[line,col]&line > col)
          cat (sprintf ("\t\\path [bidirectional] (v%d) edge (v%d);\n", line, col),
               sprintf ("\t\\path [bidirectional] (v%d) edge (v%d);\n", col, line)) ##edges drawing
        else if (!bidirectional[line,col]&adj[line,col]) 
          cat (sprintf ("\t\\path [unidirectional] (v%d) edge (v%d);\n", line, col)) ##edges drawing
      }
    }
  
  cat("\\end{tikzpicture}\n")
}