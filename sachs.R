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

datalist <- lapply(
  completepaths,
  FUN = function(fn) {
    tmp <- read.csv(fn)
    colnames(tmp) <- tolower(colnames(tmp))
    return(tmp)
  }
)

dataall <- datalist[[1]]
for (i in 2:length(datalist)) {
  dataall <- rbind(dataall, datalist[[i]])
}

Ns <- sapply(datalist, nrow)
N <- nrow(dataall)
p <- ncol(dataall)
idxs <- rep(0, length(Ns))
for (i in 2:length(Ns)) {
  idxs[i] <- Ns[i - 1] + idxs[i - 1]
}

############################  ALL DATASET

tmp <- lapply(1:10, {
  idxT <- sample(1:N, size = N / 2) #TEST
  SigmaTrain <- cor(dataall[-idxT, ])
  SigmaTest <-  cor(dataall[idxT, ])
  resultspath <- llBpath(SigmaTrain,
                         lambdas = seq(0, 2, length.out = 100),
                         eps = 1e-6)
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
  tmp <- sapply(resultspath, function(res) {
    mll(solve(res$Sigma), SigmaTest)
  })
  bidx <- which.min(tmp)
  print(bidx)
  return(resultspath[[bidx]])
})

############################### DIFFERENTE CONDITIONS
results <- lapply(1:length(idxs), function(i) {
  idx <- idxs[i] + 1:Ns[i]
  idxTest <- sample(idx, size = length(idx) / 2) #TEST
  idxTrain <- idx[ !(idx %in% idxTest) ]
  SigmaTrain <- cor(dataall[idxTrain, ])
  SigmaTest <-  cor(dataall[idxTest, ])
  resultspath <- llBpath(SigmaTrain,
                         lambdas = seq(0, 2, length.out = 100),
                         eps = 1e-6)
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
  tmp <- sapply(resultspath, function(res) {
    mll(solve(res$Sigma), SigmaTest)
  })
  bidx <- which.min(tmp)
  print(bidx)
  return(resultspath[[bidx]])
})

B <- matrix(nrow = p, ncol = p, 
            data = rowMeans(sapply(results, function(res) sign(abs(res$B)))))


############################### plotting 
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
adjmat <- sign(abs(t(B)))
colnames(adjmat) <- rownames(adjmat) <- colnames(dataall)
estgraph <- graph_from_adjacency_matrix(adjmat, diag = FALSE)
plot(estgraph,
     edge.arrow.size = 0.3,
     layout = layout_in_circle(estgraph, order = order))
igraph.to.tikz(estgraph, layout_in_circle(estgraph, order = order))
