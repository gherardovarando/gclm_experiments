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

dataall <- datalist[[1]]
for (i in 2:length(datalist)){
  dataall <- rbind(dataall, datalist[[i]])
}

########################## ESTIMATE GRAPHS FOR EACH CONDITION
results_A <- sapply(1:length(datalist), function(i) {
  D <- datalist[[i]]
  idx <- (1:length(datalist))[-i]
  DD <- datalist[[idx[1]]]
  for (j in 2:length(idx)){
    DD <- rbind(DD, datalist[[j]])
  }
  message("estimating for dataset ",i, "/", length(datalist))
  SigmaTrain <- cor(D)
  SigmaTest <-  cor(DD)
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
  
  return(B)
})
message("done")


B_A <- matrix(nrow = p, ncol = p, 
            data = rowMeans(sign(abs(results_A))))


########################################### REPETITIONS
nreps <- 100
results_B <- sapply(1:nreps, function(rep){
  idx <- sample(1:nrow(dataall), size = nrow(dataall) / 10, replace = FALSE)
  SigmaTrain <- cor(dataall[idx,])
  SigmaTest <- cor(dataall[-idx,])
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
      eps = 1e-6,
      job = 10,
      maxIter = 5000
    )
  })
  ### compute minus loglik
  tmp <- sapply(resultspath, function(res) {
    mll(solve(res$Sigma), SigmaTest)
  })
  bidx <- which.min(tmp)
  message("replicate ", rep, "/", nreps, " chosen:", bidx)
  resultspath[[bidx]]$B
})

B_B <- matrix(nrow = 11, ncol = 11, 
              data = rowMeans(sign(abs(results_B)))) 

#####################################################
results_C <- array(data = NA, dim = c(11, 11, 9, nreps))
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
     results_C[,,i, rep] <- resultspath[[bidx]]$B
   }
 }

B_C <- apply(results_C, MARGIN = c(1,2), function(x) mean(abs(sign(x))))


save(file = "proteinsignaling/results.RData", list = c("results_A", 
                                                       "results_B",
                                                       "results_C",
                                                       "B_A", "B_B", "B_C"))

