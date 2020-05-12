source("functions/util.R")
require(igraph)

nicenames <- c("Raf", "Mek$1/2$", "PLC$_\\gamma$", "PIP2", 
               "PIP3", "Erk$1/2$", "Akt", "PKA", "PKC", "p38", "JNK" )
layout <- matrix(nrow = 11, ncol = 2, byrow = TRUE,
                 data = c(10, 6, #RAF 
                          10, 3, #MEK
                           0, 7, #PLCG
                           0, 2, #PIP2
                           2, 4, #PIP3
                          10, 0, #ERK (P44.42)
                           6, 0, #AKT (PAKTS473)
                           7, 6.5, #PKA
                           5, 8, #PKC
                           5, 5, #P38
                           3, 5  #JNK
                          )) 

nreps <- 200
############# collect results
res <- array(dim = c(11,11,9,nreps))
for (i in 1:9){
  load(paste0("proteins/results",i, ".RData"))
  res[,,i,] <- results[,,] 
}


B <- apply(apply(sign(abs(res)), c(1,2,3), mean), c(1,2), function(x) {
  mean(x > 0.85)
})

print(B)
savegraphs(B, "plot/proteins/graphs/", layout = layout, names = nicenames )

# ##### GRAPH FROM SACHS et al. 
# names <- c("Raf", "Mek", "PLC", "PIP2", 
#            "PIP3", "Erk", "Akt", "PKA", "PKC", "p38", "JNK" )
# A <- matrix(nrow = 11, ncol = 11, dimnames = list(names, names),0)
# A["PIP3", "PIP2"] <- A["PIP2", "PIP3"] <- 1
# A["PLC", "PIP2"] <- 1
# A["PIP3", "PLC"] <- 1
# A['PLC', "PKC"] <- 1
# A['PIP2', "PKC"] <- 1
# A["PKC", "JNK"] <- 1
# A["PKC", "p38"] <- 1
# A["PKC", "Raf"] <- 1
# A["PIP3", "Akt"] <- 1
# A["PKA", "Akt"] <- 1
# A["PKA", "Raf"] <- 1
# A["PKA", "p38"] <- 1
# A["PKA", "Erk"] <- 1
# A["PKA", "Mek"] <- 1
# A["PKA", "JNK"] <- 1
# A["Mek", "Erk"] <- 1
# A["Raf", "Mek"] <- 1



