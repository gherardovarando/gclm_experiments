dir.create("sachs/graphs", 
           showWarnings = FALSE, recursive = TRUE)

load("sachs/results.RData")
nicenames <- c("Raf", "Mek$1/2$", "PLC$_\\gamma$", "PIP2", 
               "PIP3", "Erk$1/2$", "Akt", "PKA", "PKC", "p38", "JNK" )
layout <- matrix(nrow = 11, ncol = 2, byrow = TRUE,
                 data = c(10, 4, #RAF 
                          10, 2, #MEK
                           0, 7, #PLCG
                           0, 2, #PIP2
                           3, 6, #PIP3
                          10, 0, #ERK (P44.42)
                           5, 6, #AKT (PAKTS473)
                           6, 4, #PKA
                           5, 8, #PKC
                           6, 0, #P38
                           3, 0  #JNK
                          )) 

dir.create("sachs/graphs/A/", 
           showWarnings = FALSE, recursive = TRUE)
savegraphs(B_A, "sachs/graphs/A/", layout = layout, names = nicenames )
dir.create("sachs/graphs/B/", 
           showWarnings = FALSE, recursive = TRUE)
savegraphs(B_B, "sachs/graphs/B/", layout = layout, names = nicenames )
dir.create("sachs/graphs/C/", 
           showWarnings = FALSE, recursive = TRUE)
savegraphs(B_C, "sachs/graphs/C/", layout = layout, names = nicenames )
