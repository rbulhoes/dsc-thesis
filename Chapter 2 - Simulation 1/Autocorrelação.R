##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim1")
getwd()

##### Some packages
library(coda)

##### Analysis
e_D_TamanhoCem_AF <- read.table("TamanhoCem_AF/e_D.txt",
                                head = FALSE,
                                sep = ",")
colnames(e_D_TamanhoCem_AF) <- c("Iteration", "Site",
                                 "SS - Coord. 1", "SS - Coord. 2")

e_D_TamanhoCem_MH <- read.table("TamanhoCem_MH/e_D.txt",
                                head = FALSE,
                                sep = ",")
colnames(e_D_TamanhoCem_MH) <- c("Iteration", "Site",
                                 "MH - Coord. 1", "MH - Coord. 2")

## Autocorrelation - A4 Portrait
par(mfrow = c(2,2))
autoc_coord1_AF <- as.mcmc(x = subset(e_D_TamanhoCem_AF,
                           subset = Site == 3,
                           select = c("SS - Coord. 1")))
autocorr(x = autoc_coord1_AF)
autocorr.plot(x = autoc_coord1_AF, auto.layout = FALSE)

autoc_coord1_MH <- as.mcmc(x = subset(e_D_TamanhoCem_MH,
                                      subset = Site == 3,
                                      select = c("MH - Coord. 1")))
autocorr(x = autoc_coord1_MH)
autocorr.plot(x = autoc_coord1_MH, auto.layout = FALSE)

autoc_coord2_AF <- as.mcmc(x = subset(e_D_TamanhoCem_AF,
                           subset = Site == 3,
                           select = c("SS - Coord. 2")))
autocorr(x = autoc_coord2_AF)
autocorr.plot(x = autoc_coord2_AF, auto.layout = FALSE)

autoc_coord2_MH <- as.mcmc(x = subset(e_D_TamanhoCem_MH,
                                      subset = Site == 3,
                                      select = c("MH - Coord. 2")))
autocorr(x = autoc_coord2_MH)
autocorr.plot(x = autoc_coord2_MH, auto.layout = FALSE)

dev.off()

