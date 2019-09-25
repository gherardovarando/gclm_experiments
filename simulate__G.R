library(clggm)
library(igraph)
library(ggplot2)
source("functions/util.R")

p <- 20
B <- matrix(nrow = p, ncol = p,
            data = c(0,0,0,0,0,0,0,0,0,0,0))