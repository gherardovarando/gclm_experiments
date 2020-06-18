source("functions/util.R")

rep <- 100
N <- 1000
ks <- c(1,2,3,4)
algsel <- c("loglik", "frobenius", "pnll", "lasso", "glasso", "covthr") 

alldf <- data.frame()
itrdf <- data.frame()
for (p in 10*c(1, 2, 3, 4, 5, 10)){
     load(paste0("simulations/results_p",p,"_N", N , ".RData"))
     avgrestable <- apply(restable, MARGIN = c(1,2,5), mean, na.rm = TRUE)
     df <- expand.grid(dimnames(avgrestable))
     df <- cbind(df,Y = (apply(df, 1, function(x){
           avgrestable[x[1], x[2], x[3]]
     })))
     df1 <- df[df$stats %in% c("elapsed") & df$N == "1000",]
     df2 <- df[df$stats %in% c("iter") & df$N == "1000",]
#     df1 <- df1[df1$Algorithm %in% algsel, ]
     df1$N <- as.numeric(levels(df1$N))[as.numeric(df1$N)]

     df2$N <- as.numeric(levels(df1$N))[as.numeric(df1$N)]
     df1$p <- p
     df2$p <- p

     alldf <- rbind(alldf,df1)
     itrdf <- rbind(itrdf, df2)
}
levels(alldf$Algorithm)[levels(df$Algorithm)== "loglik"] <- "mloglik-inf" 
levels(alldf$Algorithm)[levels(df$Algorithm)== "pnll"] <- "mloglik-0.01" 
levels(alldf$Algorithm)[levels(df$Algorithm)== "frobenius"] <- "frob-inf" 

levels(itrdf$Algorithm)[levels(df$Algorithm)== "loglik"] <- "mloglik-inf" 
levels(itrdf$Algorithm)[levels(df$Algorithm)== "pnll"] <- "mloglik-0.01" 
levels(itrdf$Algorithm)[levels(df$Algorithm)== "frobenius"] <- "frob-inf" 



plotpath <- paste0("plot/simulations" )

ggplot(alldf, aes(x = p, y = Y, group = Algorithm, 
                color = Algorithm)) + 
#  facet_grid(cols = vars(k), 
#             rows = vars(stats), 
#             scales = "free_y") + 
  geom_path() + 
#  geom_hline(data = df2, 
#             aes(yintercept = Y), linetype = "dotted") +
  theme_bw() + 
#  scale_x_log10() +
#  scale_y_continuous(limits = function(x){
#    x <- x + c(-0.0005, 0.0005) * x
#    x <- c(floor(x[1] * 10),
#           ceiling(x[2] * 10)) / 10
#    return(x)
#  }, breaks = function(x) {
#    seq(x[1], x[2], length.out = 11)[c(2,4,6,8,10)]
#  }, expand = expand_scale(0,0)) + 
  theme(legend.position = "none", 
    #    axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30)) + 
  scale_y_log10(breaks = c(0.01, 1, 100)) +  
  ylab("seconds") +  
  ggsave(paste0("times", ".pdf"), path = plotpath,
         width = 3, height = 2, units = "in")




ggplot(itrdf, aes(x = p, y = Y, group = Algorithm, 
                color = Algorithm)) + 
#  facet_grid(cols = vars(k), 
#             rows = vars(stats), 
#             scales = "free_y") + 
  geom_path() + 
#  geom_hline(data = df2, 
#             aes(yintercept = Y), linetype = "dotted") +
  theme_bw() + 
#  scale_x_log10() +
#  scale_y_continuous(limits = function(x){
#    x <- x + c(-0.0005, 0.0005) * x
#    x <- c(floor(x[1] * 10),
#           ceiling(x[2] * 10)) / 10
#    return(x)
#  }, breaks = function(x) {
#    seq(x[1], x[2], length.out = 11)[c(2,4,6,8,10)]
#  }, expand = expand_scale(0,0)) + 
  theme(legend.position = "none", 
    #    axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30)) + 
  scale_y_log10() +  
  ylab("#iter") +  
  ggsave(paste0("iter", ".pdf"), path = plotpath,
         width = 3, height = 2, units = "in")


