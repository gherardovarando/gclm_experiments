args = commandArgs(trailingOnly=TRUE)  
library(ggplot2)

rep <- 100
N <- 1000
ks <- c(1,2,3,4)
ps <- 10 * (1:5)
algs <- c("loglik", "covthr")
bpath <- "simulations/"
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
  ixN    <- which(args %in% "N")  
  if (length(ixN) == 1){
     if (ixN < la){
         N <- as.numeric(args[ixN + 1]) 
         message("N set to ", N)
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
p = ps[1]
P = p
load(paste0(bpath, "results_p",p,".RData"))

avgrestable <- apply(restable[,,,paste0(p),paste0(N),], MARGIN = 1:3, mean)

df <- expand.grid(dimnames(avgrestable))
df <- cbind(df,Y = (apply(df, 1, function(x){
  avgrestable[x[1], x[2], x[3]]
})))
df$d <- as.numeric(df$k) / p
df$p <- as.numeric(p)

for (p in ps[-1]){
    
load(paste0(bpath, "results_p",p,".RData"))


avgrestable <- apply(restable[,,,paste0(p),paste0(N),], MARGIN = 1:3, mean)

dft <- expand.grid(dimnames(avgrestable))
dft <- cbind(dft,Y = (apply(df, 1, function(x){
  avgrestable[x[1], x[2], x[3]]
})))
dft$d <- as.numeric(dft$k) / p
dft$p <- as.numeric(p)
df <- rbind(df, dft) 
}

plotpath <- paste0("plot/", bpath, "N", N, "/" )
dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)


selected <- c("maxf1", "auroc", "aupr")

 
df1 <- df[df$stats %in% selected,]
df1$k <- paste0("k=",df1$k)
df2 <- data.frame(k = paste0("k=",ks), 
                  stats = "maxbacc", Y =  1 - ks / p)
ggplot(df1, aes(x = p, y = Y, group = Algorithm, 
                color = Algorithm)) + 
  facet_grid(cols = vars(k), 
             rows = vars(stats), 
             scales = "free_y") + 
  geom_path() + 
#  geom_hline(data = df2, 
#             aes(yintercept = Y), linetype = "dotted") +
  #scale_x_log10() +
  theme_bw() + 
 # ggtitle(title) + 
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
  ggsave(paste0("p", p,"N",N, ".pdf"), path = plotpath,
         width = 6, height = 3, units = "in")

  message("saved plot")  

