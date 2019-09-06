library(ggplot2)

p <- 20
d <- 1 / p
rep <- 100
dir.create("plot", showWarnings = FALSE)
dir.create("plot/simulated_ou", showWarnings = FALSE)

plotpath <- "plot/simulated_ou/"
restable <- data.frame(N = c(50, 100))
for (i in 1:nrow(restable)){
  N <- restable[i, "N"]
  dirname <- paste0("p",p,"N",N,"/", "d", d, "/")
  dir.create(paste0(plotpath, dirname),showWarnings =  FALSE, recursive = TRUE)
  path <- paste0("simulated_OU/", dirname )
  expres <- lapply(1:rep, function(ii){
    name <- paste0("S_OU_", ii, "_",
                   "p",p,"N",N, "d", d,
                   ".RData")
    load(paste0(path, name))
    #if (any(range(roc[,1]) != c(0,1) | 
    #        range(roc[,2]) != c(0,1))){
    #  warning("CHECK ROC LIMITS")
    #}
    return(list(roc = roc, auroc = auroc))
  })
  
  df <- as.data.frame(expres[[1]]$roc)
  df$g <- 1
  for (j in 2:length(expres)){

    dtmp <- as.data.frame(expres[[j]]$roc)
    dtmp$g <- j
    df <- rbind(df, dtmp)
  }
  
  
  ggplot(data=df, aes(x= FPR, y= TPR, group = g)) +
    geom_path()+
    geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") + 
    coord_fixed() + 
    ggtitle(paste("ROC"), subtitle = paste0("p=", p, ", N=", N, ", d=", d))+
    ggsave(paste0(plotpath,dirname,"roc_envelope.pdf"))
  
  restable[i, "avg_auroc"] <- mean(sapply(expres, function(x) x$auroc))

}

ggplot(data=restable, aes(x= N, y= avg_auroc)) +
  geom_path() + ylab("average auroc") + 
  ggsave(paste0(plotpath,"p",p,"d", d,"avg_roc.pdf"))

