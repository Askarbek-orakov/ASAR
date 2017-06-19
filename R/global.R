#colors of the heatmap
#names: titles, panel names
# default mgm, sp, fl, pid, txl,
# 
# library(shiny)
# library(pander)
# library(knitr)
# library(ggplot2)
# library(mmnet)
# library(RCurl)
# library(metagenomeSeq)
# library(gplots)
# library(xtable)
# library(data.table)
# library(RJSONIO)
# library(plyr)
# library(pathview)
# library(stringr)
# library(biomformat)
# library(KEGGREST)
# library(devtools)


maintitle <- "METAGENOMIC ANALYSIS by ASAR"
sel <- c(colnames(d.bm[,4]), colnames(d.bm[,5])) 
sell <- c(colnames(d.bm[,4]))
selle <- c(d.bm$usp[(nrow(d.bm)/2)])
sellec <- as.character(funtaxall$genus[(nrow(funtaxall)/2)])
sellect <- as.character(funtaxall$FUN4[(nrow(funtaxall)/2)])