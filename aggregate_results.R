args = commandArgs(trailingOnly=TRUE)  
source("functions/util.R")

rep <- 100 
Ns <- as.character(1000)
ks <- c(1,2,3,4)
Ps <- as.character(c(10, 12, 15, 20, 25, 30))
n <- length(Ns)
p <- 10
bpath <- "simulations/"
algs <- c("loglik", "frobenius", "pnll", "lasso", 
                     "glasso", "covthr") 
if (length(args) != 0){ 
  message("arguments found, parsing...")
  la <- length(args) 
  ixpath <- which(args %in% "path")
  if (length(ixpath) == 1){
     if (ixpath < la){
         bpath <- args[ixpath + 1] 
         message("base path set to ", bpath)
     }
  }
  ixp    <- which(args %in% "p")  
  if (length(ixp) == 1){
     if (ixp < la){
         p <- as.numeric(args[ixp + 1]) 
         message("p set to ", p)
     }
  }
  ixn    <- which(args %in% "N")  
  if (length(ixn) == 1){
     if (ixn < la){
         Ns <- args[ixn + 1] 
         message("N set to ", Ns)
     }
  }
  ixP    <- which(args %in% "P") 
  if (length(ixP) == 1){
     if (ixP < la){
         if (ixP > max(ixp, ixpath)){
             Ps <- (args[(ixP + 1):(la)]) 
         }else{
             Ps <- (args[ixP + 1])
         }
         message("P(s) set to ", Ps)
     }
  }
}
restable <- array(dim = c(10, length(algs), length(ks), length(Ps), length(Ns), rep), 
                  dimnames = list(stats = c("auroc", "maxacc", "maxf1", 
                                            "maxbacc", "elapsed", "aupr",
                                            "indxEmpty", "minrecall",
                                            "maxrecall", "iter"),
                                  Algorithm = algs,
                                  k = ks,
                                  P = Ps,
                                  N = Ns,
                                  rep = 1:rep), data = NA)
 
for (P in Ps){
    for (k in ks){
      for (N in Ns){
      for (i in 1:rep){
        filepath <- paste0(bpath, "p",p, "/P", P , "/k", k,"/N", N ,"/", 
                           "rep", i, ".RData" )
        load(filepath)
        P <- paste0(P)
        N <- paste0(N)
        npar <- sum(exper$B !=0 ) - p
          evals <- lapply(results[[paste0(N)]], evaluatePathB, B = exper$B[1:p,1:p])
          algs <- names(evals)
          restable["auroc" ,algs,k,P,N,i] <- sapply(evals, function(x) 
            AUROC(x$roc))
          restable["maxacc",algs,k,P,N,i] <- sapply(evals, function(x){
            ( 1 - (min(x$confusion$err) / (p*(p - 1)) ))  
          })
          restable["maxf1",algs,k,P,N,i] <- sapply(evals, function(x) 
            max(x$confusion$f1) )
          restable["maxbacc",algs,k,P,N,i] <- sapply(evals, function(x) 
            max(x$confusion$bacc) )
          restable["elapsed", algs, k, P, N, i] <- sapply(times[[paste0(N)]], 
                                                       function(x) x[3])
          restable["indxEmpty",algs,k,P,N,i] <- sapply(evals, function(x) 
            min(which(x$confusion$npar == 0 ) ) )
          restable["aupr" ,algs,k,P,N,i] <- sapply(evals, function(x) 
            AUCPR(x$confusion))
          restable["minrecall" ,algs,k,P,N,i] <- sapply(evals, function(x) 
            min(x$confusion$recall))
          restable["maxrecall" ,algs,k,P,N,i] <- sapply(evals, function(x) 
            max(x$confusion$recall))
           restable["iter" ,algs,k,P,N,i] <- sapply(evals, function(x) 
sum(x$confusion$iter ))
        }
      }
    }
    
    # if (all(restable["minrecall", , , ,] == 0)){
    #   message("min recall ok")
    # }else{
    #   message("min recall problem")
    # }
    # 
    # if (all(restable["maxrecall", , , , ] == 1)){
    #   message("max recall ok")
    # }else{
    #   message("max recall problem")
    # }
    
      message("computed p=",p, " P=", P)  
}

save(file = paste0(bpath , "results_p",p,"_N", N ,".RData" ), list = c("restable"))



