library(clggm)
library(igraph)
source("functions/util.R")

dir.create("proteinsignaling/", 
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

message("data loaded correctly")

#############################################

p < -11
nreps <- 100
results<- array(data = NA, dim = c(11, 11, 9, nreps))
 for (i in 1:9){
   message("estimating for dataset ",i, "/", length(datalist))
   D <- datalist[[i]]
   for (rep in 1:nreps){
     message("rep ",rep, "/", nreps)
     ixs <- sample(1:nrow(D), size = nrow(D)/2, replace = FALSE)
     SigmaTrain <- cor(D[ixs,])
     SigmaTest <-  cor(D[-ixs,])
     ### estimate path 
     resultspath <- llBpath(SigmaTrain, 
                            lambdas = 2*exp((1-(100:1))/10) ,
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
     results[,,i, rep] <- resultspath[[bidx]]$B
   }
 }



save(file = "proteinsignaling/results.RData", list = c( "results"))

