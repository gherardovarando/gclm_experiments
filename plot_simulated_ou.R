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

load(paste0("simulations/results_p",p,".RData"))
for (P in Ps){
    plotpath <- paste0("plot/simulations/","p",p , "/P" , P, "/" )
    dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)
   
    avgrestable <- apply(restable[,,,P,,], MARGIN = 1:4, mean)
    
    df <- expand.grid(dimnames(avgrestable))
    df <- cbind(df,Y = (apply(df, 1, function(x){
      avgrestable[x[1], x[2], x[3], x[4]]
    })))
    df$d <- as.numeric(df$k) / p
    df$N <- as.numeric(Ns)[as.numeric(df$N)]
    
    selected <- c("maxacc", "maxf1", "auroc", "aupr")
    
    df1 <- df[df$stats %in% selected,]
    df1 <- df1[df1$Algorithm %in% algsel, ]
    df1$k <- paste0("k=",df1$k)
    df2 <- data.frame(k = paste0("k=",ks), 
                      stats = "maxacc", Y =  1 - ks / p)
    ggplot(df1, aes(x = N, y = Y, group = Algorithm, 
                    color = Algorithm)) + 
      facet_grid(cols = vars(k), 
                 rows = vars(stats), 
                 scales = "free_y") + 
      geom_path() + 
      geom_hline(data = df2, 
                 aes(yintercept = Y), linetype = "dotted") +
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
            axis.text.x = element_text(angle = 30),
            legend.box.spacing = unit(0.5, "lines"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + 
      ggsave(paste0("p", p,"P",P, ".pdf"), path = plotpath,
             width = 6, height = 6, units = "in")
    
      message("saved plot p=",p, " P=", P)  
}


plotpath <- paste0("plot/simulations/" )
N <- "3000"
avgrestable <- apply(restable[,,,,N,], MARGIN = 1:4, mean)
df <- expand.grid(dimnames(avgrestable))
df <- cbind(df,Y = (apply(df, 1, function(x){
  avgrestable[x[1], x[2], x[3], x[4]]
})))
df1 <- df[df$stats %in% c("aupr"),]
df1 <- df1[df1$Algorithm %in% algsel, ]
df1$k <- paste0("k=",df1$k)

ggplot(df1, aes(x = P, y = Y, group = Algorithm, 
                color = Algorithm)) + 
  facet_grid(cols = vars(k), 
             rows = vars(stats), 
             scales = "free_y") + 
  geom_path() + 
#  geom_hline(data = df2, 
#             aes(yintercept = Y), linetype = "dotted") +
  theme_bw() + 
  scale_y_continuous(limits = function(x){
    x <- x + c(-0.0005, 0.0005) * x
    x <- c(floor(x[1] * 10),
           ceiling(x[2] * 10)) / 10
    return(x)
  }, breaks = function(x) {
    seq(x[1], x[2], length.out = 11)[c(2,4,6,8,10)]
  }, expand = expand_scale(0,0)) + 
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30)) + 
        xlab("p'") + 
  ggsave(paste0("p", p,"N",N, ".pdf"), path = plotpath,
         width = 6, height = 2, units = "in")





