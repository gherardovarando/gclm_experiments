args = commandArgs(trailingOnly=TRUE)  

library(gclm)
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

bpath <- "proteins/all/"
dir.create(bpath, showWarnings = FALSE, recursive = TRUE)

###################### arguments
k <- 2
conditions <- 1:9
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


########################## LOADING DATA 

datalist <- lapply(
  conditions,
  FUN = function(cond) {
    fn <- completepaths[cond]
    tmp <- read.csv(fn)
    colnames(tmp) <- tolower(colnames(tmp))
    return(tmp)
  }
)

p <- ncol(datalist[[1]])


all <- datalist[[1]]
if (length(datalist) > 1){
for (cond in 2:length(datalist)){
  all <- rbind(all, datalist[[cond]])
}
}

message("data loaded correctly ", dim(all))

#############################################

nreps <-200
results<- array(data = NA, dim = c(11, 11, nreps))
 for (rep in 1:nreps){
     message("  rep ",rep, "/", nreps)
     ixs <- sample(1:nrow(all), size = nrow(all)/k, replace = FALSE)
     SigmaTrain <- cov(all[-ixs,])
     dd <- diag(1 / sqrt(diag(SigmaTrain)))
     SigmaTest <-  cov(all[ixs,])
     ### estimate path 
     resultspath <- gclm.path(cov2cor(SigmaTrain), 
                             B = - 0.5 * solve(cov2cor(SigmaTrain)), 
                            lambdac = 0.01,
                            lambdas = 6 * 10^seq(-4,0, length = 100) ,
                            eps = 1e-6, job = 0, maxIter = 1000)
     ### fit MLE to all path
     resultspath <- lapply(resultspath, function(res) {
       gclm(
         cov2cor(SigmaTrain),
         B = res$B,
         C = res$C,
         lambda = 0,
         lambdac = -1,
         eps = 1e-10,
         job = 10
       )
     })
#     ### compute minus loglik
     tmp <- sapply(resultspath, function(res) {
       mll(solve(res$Sigma), dd %*% SigmaTest %*% dd)
     })
     bidx <- which.min(tmp)
     results[,, rep] <- resultspath[[bidx]]$B
   }
   save(file = paste0(bpath,"results_all.RData"), list = c( "results"))




