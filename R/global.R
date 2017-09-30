maintitle <- ""

titletab1 <- "F/T"
titletab2 <- "F/M"
titletab3 <- "T/M"
titletab4 <- "KEGG Pathways Heatmap"
titletab5 <- "KEGG Pathway"
titletab6 <- "Settings"
titletab7 <- "Metadata"
titleForEmtyDataPopup <- "Empty selection!"
textForEmptyDataPopup <- "Please select some function or taxon."
titleIntfuntax <- "Heatmap cannot be generated"
textIntfuntaxInPlot1 <- "There are no data to display for this function and taxonomy in your dataset. Please try to change either function or taxonomy selections."
textIntfuntaxInPlot2 <- "There are no data to display for this function and taxonomy in your dataset. Please try to change either function or taxonomy selections."
textIntfuntaxInPlot3 <- "There are no data to display for this function and taxonomy in your dataset. Please try to change either function or taxonomy selections."
titleForDimErrorPopup <- "Wrong Heatmap dimensions!"
textForDimErrorPopup <- "Please select other function or taxon or go one level higher."
titleForNonMatchKostat <- "Error with annotation files!"
textForNonMatchKostat <- "There is inconsistency within KEGG Orthology Annotation file."
titleForSavingMetadata <- "Saving your metadata..."
textForSavingMetadata <- "Your changes will be applied after you restart the app. Saving process may take a few seconds. Please, restart the app!"
textForDownloadingMetadata <- "KEGG Map is saved in your working directory."
ColNameSelectorTitle <- "Please select a column of metadata to use as metagenome names"
metagenomeone <-"Select multiple metagenomes"
metagenome1selected <- c(14,15)
metagenometwo <- "Select metagenome"
metagenome2selected <- 14

plot1Title <-'F vs. T heatmap'
plot2Title <-'F vs. M heatmap'
plot3Title <-'T vs. M heatmap'
plot4Title <-'KEGG pathway heatmap'


taxone <- "Taxonomy Selection level"
tax1n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain", "root" = "toplevel")
tax1selected <- "genus"
taxthree <- "Select taxons"
taxtwo <- "Taxonomy Aggregation level"
tax2n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain")
tax2selected <- "usp"

funcone <- "Function Selection level"
func1n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4", "root" = "toplevel")
func1selected <- "FUN4"
functhree <- "Select function"
functwo <- "Function Aggregation level"
func2n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4")
func2selected <- "ufun"

pathwayone <- "Input Pathway ID"

#colname for sample name pickup.
#Color Palette for Heatmaps
#Settings
labelColPal <- "Select color palette for heatmaps"
currentPalette <- "Blues"
colName <- "Group"
if(file.exists("Settings.Rdata")){
  load("Settings.Rdata")
}
