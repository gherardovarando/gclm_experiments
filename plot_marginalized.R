args = commandArgs(trailingOnly=TRUE)  
library(ggplot2)

rep <- 100
N <- 1000
ks <- c(1,2,3,4)
Ps <- as.character(c(10, 12, 15, 20, 25, 30))
p <- 10
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
load(paste0(bpath, "results_p",p,"_N", N ,".RData"))

plotpath <- paste0("plot/", bpath, "p",p , "/N", N, "/" )
dir.create(plotpath, showWarnings = FALSE, recursive = TRUE)

avgrestable <- apply(restable[,,,,paste0(N),], MARGIN = c(1,2,4), mean)

df <- expand.grid(dimnames(avgrestable))
df <- cbind(df,Y = (apply(df, 1, function(x){
  avgrestable[x[1], x[2], x[3]]
})))
#df$d <- as.numeric(df$k) / p
df$P <- df$P

selected <- c("maxacc", "maxf1", "auroc", "aupr")


levels(df$Algorithm)[levels(df$Algorithm)== "loglik"] <- "mloglik-inf" 
levels(df$Algorithm)[levels(df$Algorithm)== "pnll"] <- "mloglik-0.01" 
levels(df$Algorithm)[levels(df$Algorithm)== "frobenius"] <- "frob-inf" 

df1 <- df[df$stats %in% selected,]
#df1$k <- paste0("k=",df1$k)
df2 <- data.frame(k = paste0("k=",ks), 
                  stats = "maxacc", Y =  1 - ks / p)
print(nrow(df1))
ggplot(df1, aes(x = P, y = Y, group = Algorithm, 
                color = Algorithm)) + 
  facet_grid( #cols = vars(k), 
             rows = vars(stats), 
             scales = "free_y") + 
  geom_path() + 
#  geom_hline(data = df2, 
#             aes(yintercept = Y), linetype = "dotted") +
  #scale_x_log10() +
  theme_bw() + 
 # ggtitle(title) + 
  scale_y_continuous(limits = function(x){
    x <- x + c(-0.00005, 0.00005) * x
    x <- c(floor(x[1] * 60),
      ceiling(x[2] * 60)) / 60
    return(x)
  }, breaks = function(x) {
    round(seq(x[1], x[2], length.out = 7)[c(2,4,6)],2)
  }, expand = expand_scale(0,0)) + 
  xlab("p") + 
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30),
        legend.box.spacing = unit(0.5, "lines"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + 
  ggsave(paste0("p", p,"N",N, ".pdf"), path = plotpath,
         width = 3, height = 4.5, units = "in")

  message("saved plot")  

