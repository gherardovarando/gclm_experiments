library(clggm)
require(igraph)
source("functions/util.R")

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

## set seed 
set.seed(1)
########################## ESTIMATE GRAPH
results <- lapply(1:length(datalist), function(i) {
  D <- datalist[[i]]
  message("estimating for dataset ",i, "/", length(datalist))
  replicate(100, {
    idxTest <- sample(nrow(D), size = nrow(D) / 2) #TEST
    SigmaTrain <- cor(D[-idxTest, ])
    SigmaTest <-  cor(D[idxTest, ])
    ### estimate path 
    resultspath <- llBpath(SigmaTrain,
                           lambdas = seq(0, 3, length.out = 100),
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
    return(sign(abs(resultspath[[bidx]]$B)))
  })
})
message("done")

averages <- lapply(results, function(res){
  matrix(nrow = p, ncol = p, data = apply(res, MARGIN = c(1,2), FUN = mean))
})
B <- matrix(nrow = p, ncol = p, 
            data = rowMeans(sapply(averages, function(x) x > 0.2 )))

#### saving results

save(file = "sachsresults.RData", list = c("B", "results", "averages"))

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
thrs <- c(0, sort(unique(c(B))))
for (thr in thrs){
  adjmat <- sign(abs(t(B > thr)))
  colnames(adjmat) <- rownames(adjmat) <- colnames(datalist[[1]])
  estgraph <- graph_from_adjacency_matrix(adjmat, diag = FALSE)
  plot(estgraph,
       edge.arrow.size = 0.3,
       layout = layout_in_circle(estgraph, order = order))
  
  #### saving graph to tkiz format
  message("saving graph to sachsGraph_thr**.txt")
  igraph.to.tikz(estgraph, layout_in_circle(estgraph, order = order), 
                 file = paste0("sachsGraph_thr",thr,  ".txt"))
  
}


dataall <- datalist[[1]]
for (i in 2:length(datalist)){
  dataall <- rbind(dataall, datalist[[i]])
}
stab <- llBstabilitypath(dataall, job = 11, eps = 1e-6,
                         lambdas = seq(0,2,length.out = 100))

plot.new()
plot.window(xlim = c(0,100), ylim = c(0,1))
for (i in 1:p){
  for (j in 1:p){
    lines(stab[i,j, 100: 1], col = "black")         
  }
}
