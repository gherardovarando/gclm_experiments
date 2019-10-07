args = commandArgs(trailingOnly=TRUE)  

library(clggm)
library(igraph)
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

###################### arguments
conditions <- 1:9
k <- 2
if (length(args) != 0){ 
  message("arguments found, parsing...")
  la <- length(args) 
  ixpath <- which(args %in% "path")
  if (length(ixpath) == 1){
     if (ixpath < la){
         bpath <- args[ixpath + 1] 
         message("out path set to ", bpath)
     }
  } 
  ixk <- which(args %in% "k")
  if (length(ixk) == 1){
     if (ixk < la){
         k <- as.numeric(args[ixk + 1]) 
         message("number of folds set to ", k)
     }
  }
  ixcond    <- which(args %in% "condition")  
  if (length(ixcond) == 1){
     if (ixcond < la){
         if (ixcond > max(ixpath, ixk)){
             conditions <- as.numeric(args[(ixcond + 1):(la)]) 
         }else{
             conditions <- as.numeric(args[ixcond + 1])
         }
         message("Conditions set to ", conditions)
     }
  }
}

######################### creating dir

dir.create(bpath, showWarnings = FALSE, recursive = TRUE)
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

nreps <-200
results<- array(data = NA, dim = c(11, 11, nreps))
for (i in conditions){
   message("estimating for condition ",i, "/", length(conditions))
   D <- datalist[[i]]
   for (rep in 1:nreps){
     message("condition ", i, "  rep ",rep, "/", nreps)
     ixs <- sample(1:nrow(D), size = nrow(D)/k, replace = FALSE)
     SigmaTrain <- cor(D[ixs,])
     SigmaTest <-  cor(D[-ixs,])
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
     results[,, rep] <- resultspath[[bidx]]$B
   }
   save(file = paste0(bpath,"results", i,".RData"), list = c( "results"))
}




