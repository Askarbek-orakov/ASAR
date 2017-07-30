#colors of the heatmap

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

maintitle <- ""
load("mdt.Rdata")
metagenomeone <-"Choose metagenome samples"
#metagenome1n <- c(colnames(d.bm[,-c(1:3)]))
metagenome1n <- setNames(c(colnames(d.bm[,-c(1:3)])), mdt$MGN)
metagenome1selected <- c(colnames(d.bm[,5]), colnames(d.bm[,6]))

metagenometwo <- "Choose one metagenome sample"
#metagenome2n <- c(colnames(d.bm[,-c(1:3)]))
metagenome2n <- setNames(c(colnames(d.bm[,-c(1:3)])), mdt$MGN)
metagenome2selected <- c(colnames(d.bm[,4]))

taxone <- "Choose taxlevel 1"
tax1n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain", "Toplevel" = "toplevel")
tax1selected <- "genus"

taxtwo <- "Aggregation taxlevel 2"
tax2n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain")
tax2selected <- "usp"

taxthree <- "Select taxon"
functhree <- "Select function"
pathwayone <- "Input Pathway ID"

funcone <- "Choose funLevel 1"
func1n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4", "Toplevel/All" = "toplevel")
func1selected <- "FUN4"

functwo <- "Choose funLevel 2"
func2n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4")
func2selected <- "ufun"

specieone <- "Choose Specie"
specie1n <- as.vector(unique(d.bm[,"usp"]))
specie1sellected <- c(d.bm$usp[(nrow(d.bm)/2)])

#Color Palette for Heatmaps
currentPalette <- "Blues"
#Settings
set_taxone <- "taxlevel 1"
set_taxtwo <- "taxlevel 2"
set_funcone <- "funlevel 1"
set_functwo <- "funlevel 2"
if(file.exists("Settings.Rdata")){
  load("Settings.Rdata")
}
