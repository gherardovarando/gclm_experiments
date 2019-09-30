library(ggplot2)
source("functions/util.R")

rep <- 100
Ns <- as.character(c(50, 100, 200, 300, 400, 500,
                     600, 700, 800, 900, 
                     1000, 2000, 3000, 4000, 5000))
ks <- c(1,2,3,4)
Ps <- as.character(c(10, 12, 15, 20))
n <- length(Ns)
p <- 10
restable <- array(dim = c(9, length(algs), length(ks), length(Ps), length(Ns), rep), 
                  dimnames = list(stats = c("auroc", "maxacc", "maxf1", 
                                            "maxbacc", "elapsed", "aucpr",
                                            "indxEmpty", "minrecall",
                                            "maxrecall"),
                                  Algorithm = algs,
                                  k = ks,
                                  P = Ps,
                                  N = Ns, 
                                  rep = 1:rep), data = NA)
for (P in Ps){
    algs <- c("loglik", "lasso", "glasso")
    plotpath <- paste0("plot/simulations/","p",p , "P" , P, "/" )
    dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)
    for (k in ks){
      for (i in 1:rep){
        filepath <- paste0("simulations/", "p",p, "/P", P , "/k", k,"/"
                           , "rep", i, ".RData" )
        load(filepath)
        P <- paste0(P)
        for (N in Ns){
          npar <- sum(exper$B !=0 ) - p
          evals <- lapply(results[[paste0(N)]], evaluatePathB, B = exper$B[1:p,1:p])
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
          restable["aucpr" ,algs,k,P,N,i] <- sapply(evals, function(x) 
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
    
    avgrestable <- apply(restable[,,,P,,], MARGIN = 1:4, mean)
    
    df <- expand.grid(dimnames(avgrestable))
    df <- cbind(df,Y = (apply(df, 1, function(x){
      avgrestable[x[1], x[2], x[3], x[4]]
    })))
    df$d <- as.numeric(df$k) / p
    df$N <- as.numeric(Ns)[as.numeric(df$N)]
    
    selected <- c("maxacc", "maxf1", "auroc", "aucpr")
    
    df1 <- df[df$stats %in% selected,]
    df1$k <- paste0("k=",df1$k)
    df2 <- data.frame(k = paste0("k=",ks), 
                      stats = "maxacc", Y =  1 - ks / p)
    ggplot(df1, aes(x = N, y = Y, group = Algorithm, 
                    color = Algorithm)) + 
      facet_grid(cols = vars(k), 
                 rows = vars(stats), 
                 scales = "free_y") + 
      geom_path() + 
#      geom_hline(data = df2, 
#                 aes(yintercept = Y), linetype = "dotted") +
      scale_x_log10() +
      theme_bw() + 
      scale_y_continuous(limits = function(x){
        x <- x + c(-0.0005, 0.0005) * x
        x <- c(floor(x[1] * 10),
          ceiling(x[2] * 10)) / 10
        return(x)
      }, breaks = function(x) {
        seq(x[1], x[2], length.out = 11)[c(2,4,6,8,10)]
      }, expand = expand_scale(0,0)) + 
      theme(legend.position = "bottom", 
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 30)) + 
      ggsave(paste0("p", p,"P",P, ".pdf"), path = plotpath,
             width = 6, height = 6, units = "in")
    
}

plotpath <- paste0("plot/simulations/","p",p , "/" )
N <- "1000"
avgrestable <- apply(restable[,,,,N,], MARGIN = 1:4, mean)
df <- expand.grid(dimnames(avgrestable))
df <- cbind(df,Y = (apply(df, 1, function(x){
  avgrestable[x[1], x[2], x[3], x[4]]
})))
df1 <- df[df$stats %in% selected,]
df1$k <- paste0("k=",df1$k)

ggplot(df1, aes(x = P, y = Y, group = Algorithm, 
                color = Algorithm)) + 
  facet_grid(cols = vars(k), 
             rows = vars(stats), 
             scales = "free_y") + 
  geom_path() + 
  theme_bw() + 
  scale_y_continuous(limits = function(x){
    x <- x + c(-0.0005, 0.0005) * x
    x <- c(floor(x[1] * 10),
           ceiling(x[2] * 10)) / 10
    return(x)
  }, breaks = function(x) {
    seq(x[1], x[2], length.out = 11)[c(2,4,6,8,10)]
  }, expand = expand_scale(0,0)) + 
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30)) + 
  ggsave(paste0("p", p,"N",N, ".pdf"), path = plotpath,
         width = 6, height = 6, units = "in")





