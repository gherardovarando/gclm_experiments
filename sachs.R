library(clggm)
require(igraph)

filenames <- c("cd3cd28", "cd3cd28icam2", "cd3cd28+aktinhib", 
               "cd3cd28+g0076", "cd3cd28+psitect",  "cd3cd28+u0126",
               "cd3cd28+ly", "pma", "b2camp")
extension <- ".csv"
basepath <- "data/Sachs/Data Files/"
completepaths <- paste0(basepath, 1:9, ". ", filenames, extension)

datalist <- lapply(completepaths, FUN = function(fn){
  tmp <- read.csv(fn)
  colnames(tmp) <- tolower(colnames(tmp))
  return(tmp)
})

dataall <- datalist[[1]]
for (i in 2:length(datalist)){
  dataall <- rbind(dataall, datalist[[i]])
}

p <- ncol(dataall)
Sigma <- cor(dataall)
B0 <- -0.5 * diag(p) %*% solve(Sigma)
results <- pnllbc(Sigma, B = B0, C = diag(p), C0 = diag(p), eps = 1e-15,
                  maxIter = 100, intitr = 100, lambda = max(abs(Sigma)) / 35, 
                  lambdac = 0.005, job = 0 )
results$B

adjmat <- sign(abs(t(results$B)))
colnames(adjmat) <- rownames(adjmat) <- colnames(dataall)
estgraph <- graph_from_adjacency_matrix(adjmat) 
plot(estgraph, edge.arrow.size = 0.3, layout = layout_in_circle)

resultspath <- llBpath(Sigma, eps = 1e-7, maxIter = 5000, job = 11, 
                       lambdas = seq(0,1,length.out = 10))
sapply(resultspath, function(res) sum(res$B!=0))
sapply(resultspath, function(res) res$lambda)
sapply(resultspath, function(res) res$iter)
sapply(resultspath, function(res) res$diff)

