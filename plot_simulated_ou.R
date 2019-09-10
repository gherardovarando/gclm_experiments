library(ggplot2)
source("functions/util.R")

p <- 10
k <- 2
rep <- 100
Ns <- c(50, 100, 200, 300, 400, 500, 1000)
dir.create("plot", showWarnings = FALSE)
dir.create("plot/simulated_ou", showWarnings = FALSE)

plotpath <- "plot/simulated_ou/"
restable <- data.frame(N = Ns)
for (i in 1:nrow(restable)){
  N <- restable[i, "N"]
  path <- paste0("p",p,"/k",k,"/N", N, "/")
  filepath <- paste0("simulated_OU/identityC/", path )
  dir.create(paste0(plotpath, path),showWarnings =  FALSE, recursive = TRUE)
  exps <- as.data.frame(t(sapply(1:rep, function(ii){
    name <- paste0("rep", ii, ".RData")
    load(paste0(filepath, name))
    evals <- lapply(results, evaluatePathB, B = exper$B)
    AUROCS <- abs(sapply(evals, function(x) AUROC(x$roc)))
    names(AUROCS) <- paste0("AUROC.", names(evals)) 
    MINERR <- sapply(evals, function(x) min(x$confusion$err))
    names(MINERR) <- paste0("MINERR.", names(evals)) 
    NPAR <- c(NPAR = sum(exper$B != 0) - p)
    return(c(AUROCS, MINERR, NPAR))
    ### check ROC
    })))
  
    restable$AUROC.loglikL1B[i] <- mean(exps$AUROC.loglikL1B)
    restable$AUROC.lsL1B[i] <- mean(exps$AUROC.lsL1B)
    restable$AUROC.lasso[i] <- mean(exps$AUROC.lasso)
    restable$MINERR.loglikL1B[i] <- mean(exps$MINERR.loglikL1B / exps$NPAR)
    restable$MINERR.lsL1B[i] <- mean(exps$MINERR.lsL1B / exps$NPAR)
    restable$MINERR.lasso[i] <- mean(exps$MINERR.lasso / exps$NPAR)
}

df <- data.frame(AUROC = c(restable$AUROC.loglikL1B, 
                           restable$AUROC.lsL1B,
                           restable$AUROC.lasso), 
                 algorithm = c(rep("loglikL1B", 7),
                               rep("lsL1B", 7),
                               rep("lasso", 7)),
                 N = rep(Ns, 3))

ggplot(df, aes(x = N, y = AUROC, group = algorithm, color = algorithm)) + 
  geom_path()



df <- data.frame(MINERR = c(restable$MINERR.loglikL1B, 
                           restable$MINERR.lsL1B,
                           restable$MINERR.lasso), 
                 algorithm = c(rep("loglikL1B", 7),
                               rep("lsL1B", 7),
                               rep("lasso", 7)),
                 N = rep(Ns, 3))

ggplot(df, aes(x = N, y = MINERR, group = algorithm, color = algorithm)) + 
  geom_path()
