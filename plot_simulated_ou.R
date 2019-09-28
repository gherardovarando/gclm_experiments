library(ggplot2)
source("functions/util.R")

rep <- 100
Ns <- as.character(c(50, 100, 200, 300, 400, 500,
                     600, 700, 800, 900, 
                     1000, 2000, 3000, 4000, 5000))
ks <- c(1,2,3,4)
n <- length(Ns)
p <- 10
for (typeC in c("knowC")){
  for (typeB in c("general")){
    algs <- c("loglik", "frobenius", "lasso", "glasso")
    restable <- array(dim = c(9, length(algs), length(ks), length(Ns), rep), 
                      dimnames = list(stats = c("auroc", "maxacc", "maxf1", 
                                                "maxbacc", "elapsed", "aucpr",
                                                "indxEmpty", "minrecall",
                                                "maxrecall"),
                                      algorithm = algs,
                                      k = ks,
                                      N = Ns, 
                                      rep = 1:rep), data = NA)
    
    plotpath <- paste0("plot/simulations/", typeC, "/", typeB, "/p",p,"/")
    dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)
    for (k in ks){
      for (i in 1:rep){
        filepath <- paste0("simulations/", typeC, "/", typeB, "/p",p,"/k",k,"/"
                           , "rep", i, ".RData" )
        load(filepath)
        for (N in Ns){
          npar <- sum(exper$B !=0 ) - p
          evals <- lapply(results[[paste0(N)]], evaluatePathB, B = exper$B)
          restable["auroc" ,algs,k,N,i] <- sapply(evals, function(x) 
            AUROC(x$roc))

          restable["maxacc",algs,k,N,i] <- sapply(evals, function(x){
            ( 1 - (min(x$confusion$err) / (p*(p - 1)) ))  
          })
          restable["maxf1",algs,k,N,i] <- sapply(evals, function(x) 
            max(x$confusion$f1) )
          restable["maxbacc",algs,k,N,i] <- sapply(evals, function(x) 
            max(x$confusion$bacc) )
          restable["elapsed", algs, k, N, i] <- sapply(times[[paste0(N)]], 
                                                       function(x) x[3])
          restable["indxEmpty",algs,k,N,i] <- sapply(evals, function(x) 
            min(which(x$confusion$npar == 0 ) ) )
          restable["aucpr" ,algs,k,N,i] <- sapply(evals, function(x) 
            AUCPR(x$confusion))
          restable["minrecall" ,algs,k,N,i] <- sapply(evals, function(x) 
            min(x$confusion$recall))
          restable["maxrecall" ,algs,k,N,i] <- sapply(evals, function(x) 
            max(x$confusion$recall))
        }
      }
    }
    
    if (all(restable["minrecall", , , ,] == 0)){
      message("min recall ok")
    }else{
      message("min recall problem")
    }
    
    if (all(restable["maxrecall", , , , ] == 1)){
      message("max recall ok")
    }else{
      message("max recall problem")
    }
    
    avgrestable <- apply(restable, MARGIN = 1:4, mean)
    
    df <- expand.grid(dimnames(avgrestable))
    df <- cbind(df,Y = (apply(df, 1, function(x){
      avgrestable[x[1], x[2], x[3], x[4]]
    })))
    df$d <- as.numeric(df$k) / p
    df$N <- as.numeric(Ns)[as.numeric(df$N)]
    df$k <- paste0("k=",df$k)
    
    selected <- c("maxacc", "maxf1", "auroc", "aucpr")
    
    df1 <- df[df$stats %in% selected,]
    ggplot(df1, aes(x = N, y = Y, group = algorithm, 
                    color = algorithm)) + 
      facet_grid(cols = vars(k), rows = vars(stats), scales = "free_y") + 
      geom_path() + 
      scale_x_log10() +
      theme(legend.position = "bottom", axis.title.y = element_blank()) + 
      ggsave(paste0(typeC, typeB, "p", p, ".pdf"), path = plotpath,
             width = 6, height = 6, units = "in")
    
  }
}




