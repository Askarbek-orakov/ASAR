maintitle <- ""

titleForDimErrorPopup <- "Wrong Heatmap dimensions!"
textForDimErrorPopup <- "Please select other function or taxon or go one level higher."
metagenomeone <-"Choose metagenome samples"
metagenome1selected <- c(14,15)
metagenometwo <- "Choose one metagenome sample"
metagenome2selected <- 14

taxone <- "Choose taxlevel 1"
tax1n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain", "Toplevel" = "toplevel")
tax1selected <- "genus"
taxthree <- "Select taxon"
taxtwo <- "Aggregation taxlevel 2"
tax2n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain")
tax2selected <- "usp"

funcone <- "Choose funLevel 1"
func1n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4", "Toplevel/All" = "toplevel")
func1selected <- "FUN4"
functhree <- "Select function"
functwo <- "Choose funLevel 2"
func2n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4")
func2selected <- "ufun"

pathwayone <- "Input Pathway ID"

#colname for sample name pickup.
#Color Palette for Heatmaps
#Settings
set_taxone <- "taxlevel 1"
set_taxtwo <- "taxlevel 2"
set_funcone <- "funlevel 1"
set_functwo <- "funlevel 2"
currentPalette <- "Blues"
colName <- "Group"
if(file.exists("Settings.Rdata")){
  load("Settings.Rdata")
}
