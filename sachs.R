library(clggm)
library(igraph)
source("functions/util.R")

dir.create("resultsSachs/graphs", 
           showWarnings = FALSE, recursive = TRUE)
filenames <- c(
  "cd3cd28",
  "cd3cd28icam2",
  "cd3cd28+aktinhib",
  "cd3cd28+g0076",
  "cd3cd28+psitect",
  "cd3cd28+u0126",
  "cd3cd28+ly",
  "pma",
  "b2camp"
)
extension <- ".csv"
basepath <- "data/Sachs/Data Files/"
completepaths <- paste0(basepath, 1:9, ". ", filenames, extension)
order <- c(
  "pka",
  "pjnk",
  "p38",
  "pip3",
  "praf",
  "pakts473",
  "p44.42",
  "pkc",
  "pip2",
  "plcg",
  "pmek"
)
########################## LOADING DATA 

datalist <- lapply(
  completepaths,
  FUN = function(fn) {
    tmp <- read.csv(fn)
    colnames(tmp) <- tolower(colnames(tmp))
    return(tmp)
  }
)

p <- ncol(datalist[[1]])

dataall <- datalist[[1]]
for (i in 2:length(datalist)){
  dataall <- rbind(dataall, datalist[[i]])
}

########################## ESTIMATE GRAPHS
results <- lapply(1:length(datalist), function(i) {
  D <- datalist[[i]]
  message("estimating for dataset ",i, "/", length(datalist))
  SigmaTrain <- cor(D)
  SigmaTest <-  cor(dataall)
  ### estimate path 
  resultspath <- llBpath(SigmaTrain,
                           lambdas = 2*exp((1-(50:1))/5) ,
                           eps = 1e-6, job = 11, maxIter = 5000)
  ### fit MLE to all path
  resultspath <- lapply(resultspath, function(res) {
      proxgradllB(
        SigmaTrain,
        B = res$B,
        C = res$C,
        lambda = 0,
        eps = 1e-10,
        job = 10
      )
    })
  ### compute minus loglik
  tmp <- sapply(resultspath, function(res) {
      mll(solve(res$Sigma), SigmaTest)
    })
  bidx <- which.min(tmp)
  B <- resultspath[[bidx]]$B
  adjmat <- sign(abs(t(B)))
  colnames(adjmat) <- rownames(adjmat) <- colnames(datalist[[1]])
  estgraph <- graph_from_adjacency_matrix(adjmat, diag = FALSE)
  plot(estgraph,
         edge.arrow.size = 0.3,
         layout = layout_in_circle(estgraph, order = order))
    
  #### saving graph to tkiz format
  message("saving graph to graphs folder")
  igraph.to.tikz(estgraph, layout_in_circle(estgraph, order = order), 
                   file = paste0("resultsSachs/graphs/sachsGraph_condition", i,  ".txt"))#
  return(B)
})
message("done")


B <- matrix(nrow = p, ncol = p, 
            data = rowMeans(sapply(results, function(x) sign(abs(x))  )))

#### saving results

save(file = "resultsSachs/sachsresults.RData", list = c("B", "results", "averages"))

############################### plotting and saving
order <- c(
  "pka",
  "pjnk",
  "p38",
  "pip3",
  "praf",
  "pakts473",
  "p44.42",
  "pkc",
  "pip2",
  "plcg",
  "pmek"
)
thrs <- c(sort(unique(c(B))))
for (thr in thrs){
  adjmat <- sign(abs(t(B > thr)))
  colnames(adjmat) <- rownames(adjmat) <- colnames(datalist[[1]])
  estgraph <- graph_from_adjacency_matrix(adjmat, diag = FALSE)
  plot(estgraph,
       edge.arrow.size = 0.3,
       layout = layout_in_circle(estgraph, order = order))
  
  #### saving graph to tkiz format
  message("saving graph to graphs folder")
  igraph.to.tikz(estgraph, layout_in_circle(estgraph, order = order), 
                 file = paste0("resultsSachs/graphs/sachsGraph_thr",thr,  ".txt"))
  
}

