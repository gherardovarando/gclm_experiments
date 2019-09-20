library(ggplot2)
source("functions/util.R")

rep <- 20
Ns <- as.character(c(50, 100, 200, 300, 400, 500,
                     600, 700, 800, 900, 
                     1000))
ks <- c(1)
n <- length(Ns)
p <- 20
algs <- c("loglik", "frobenius", "lasso", "glasso")
restable <- array(dim = c(5, length(algs), length(ks), length(Ns), rep), 
                  dimnames = list(stats = c("AUROC", "MINERR", "MAXF1", 
                                            "MAXBACC", "ELAPSED"),
                                  algorithm = algs,
                                  k = ks,
                                  N = Ns, 
                                  rep = 1:rep), data = NA)

plotpath <- paste0("plot/simulated_ou/randomC/", "p",p,"/")
dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)
for (k in ks){
  for (i in 1:rep){
    filepath <- paste0("simulated_OU/randomC/", "p",p,"/k",k,"/"
                       , "rep", i, ".RData" )
    load(filepath)
    for (N in Ns){
      npar <- sum(exper$B !=0 ) - p
      evals <- lapply(results[[paste0(N)]], evaluatePathB, B = exper$B)
      restable["AUROC" ,algs,k,N,i] <- sapply(evals, function(x) 
                                                        AUROC(x$roc))
      restable["MAXACC",algs,k,N,i] <- sapply(evals, function(x) 
                                                   (1 - min(x$confusion$err) / npar))
      restable["MAXF1",algs,k,N,i] <- sapply(evals, function(x) 
                                                    max(x$confusion$f1) )
      restable["MAXBACC",algs,k,N,i] <- sapply(evals, function(x) 
        max(x$confusion$bacc) )
      
      restable["ELAPSED", algs, k, N, i] <- sapply(times, function(x) x[3])
    }
  }
}


avgrestable <- apply(restable, MARGIN = 1:4, mean)

df <- expand.grid(dimnames(avgrestable)[-1])
df <- cbind(df,t(apply(df, 1, function(x){
  avgrestable[,x[1], x[2], x[3]]
})))
df$d <- as.numeric(df$k) / p
df$N <- as.numeric(Ns)[as.numeric(df$N)]


ggplot(df, aes(x = N, y = AUROC, group = paste0(algorithm,k), color = algorithm, 
               linetype = k )) + 
  geom_path() + 
  ggsave("AUROC.pdf", path = plotpath,
         width = 10, height = 10, units = "cm")


ggplot(df, aes(x = N, y = MAXACC, group = paste0(algorithm,k), color = algorithm, 
               linetype = k )) + 
  geom_path() + ggsave("MAXACC.pdf", path = plotpath,
                       width = 10, height = 10, units = "cm")

ggplot(df, aes(x = N, y = MAXF1, group = paste0(algorithm,k), color = algorithm, 
               linetype = k )) + 
  geom_path() + ggsave("MAXF1.pdf", path = plotpath,
                       width = 10, height = 10, units = "cm")

ggplot(df, aes(x = N, y = MAXBACC, group = paste0(algorithm,k), color = algorithm, 
               linetype = k )) + 
  geom_path() + ggsave("MAXBACC.pdf", path = plotpath,
                       width = 10, height = 10, units = "cm")





