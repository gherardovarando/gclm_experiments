source("functions/util.R")

rep <- 100
N <- 1000
ks <- c(1,2,3,4)
P <- 10
p <- 10
algs <- c("loglik", "frobenius", "lasso",  
                     "glasso", "covthr")
restable <- array(dim = c(100, length(algs), length(ks), rep), 
                  dimnames = list(recall = 1:100,
                                  Algorithm = algs,
                                  k = ks,
                                  rep = 1:rep), data = NA)
 
    plotpath <- paste0("plot/simulations/","p",p , "/P" , P, "/N",N,"/"  )
    dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)
    for (k in ks){
      for (i in 1:rep){
        filepath <- paste0("simulations/", "p",p, "/P", P , "/k", k,"/N", N, "/", 
                       "rep", i, ".RData" )
        load(filepath)
        P <- paste0(P)
        evals <- lapply(results[[paste0(N)]], evaluatePathB, B = exper$B[1:p,1:p])
        algs <- names(evals)
        for (a in algs){
            restable[1:100, a, k,  i] <- evals[[a]]$cpr$y
        }
      }
     }   
    
    avgrestable <- apply(restable, MARGIN = 1:3, mean, na.rm = TRUE)
    nna <- apply(restable, MARGIN = c(2,3), function(x) sum(is.na(x)))
    print(nna) 
    df <- expand.grid(dimnames(avgrestable))
    df <- cbind(df,precision = (apply(df, 1, function(x){
      avgrestable[x[1], x[2], x[3]]
    })))
    df <- df[df$Algorithm %in% algs, ]
    df$k <- paste0("k=",df$k)
    df$recall <- seq(0,1,length.out = 100)
    ggplot(df, aes(x = recall, y = precision, group = Algorithm, 
                    color = Algorithm)) + 
      facet_grid(cols = vars(k), 
                 scales = "free_y") + 
      geom_path() + 
      theme_bw() + 
      theme(legend.position = "none", 
            axis.text.x = element_text(angle = 30),
            legend.box.spacing = unit(0.5, "lines"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + 
      ggsave(paste0("averageCPR.pdf"), path = plotpath,
             width = 6, height = 2, units = "in")
    
      message("saved plot average CPR")  



