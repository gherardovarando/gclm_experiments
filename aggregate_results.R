library(ggplot2)
source("functions/util.R")

rep <- 100
Ns <- as.character(c(50, 100, 200, 300, 400, 500,
                     600, 700, 800, 900, 
                     1000, 2000, 3000, 4000, 5000))
ks <- c(1,2,3,4)
Ps <- as.character(c(10, 12, 15, 20, 25, 30, 35, 40))
n <- length(Ns)
p <- 10
algs <- c("loglik", "frobenius", "lasso", "lassoc", 
                     "glasso", "covthr")
algsel <- c("loglik", "frobenius", "lasso", "glasso", "covthr") 
restable <- array(dim = c(9, length(algs), length(ks), length(Ps), length(Ns), rep), 
                  dimnames = list(stats = c("auroc", "maxacc", "maxf1", 
                                            "maxbacc", "elapsed", "aupr",
                                            "indxEmpty", "minrecall",
                                            "maxrecall"),
                                  Algorithm = algs,
                                  k = ks,
                                  P = Ps,
                                  N = Ns, 
                                  rep = 1:rep), data = NA)
 
for (P in Ps){
    for (k in ks){
      for (i in 1:rep){
        filepath <- paste0("simulations/", "p",p, "/P", P , "/k", k,"/"
                           , "rep", i, ".RData" )
        load(filepath)
        P <- paste0(P)
        for (N in Ns){
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

save(paste0(file = "simulations/results_p",p,".RData" ), list = c("restable"))



