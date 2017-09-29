_Askarbek Orakov, Nazgul Sakenova, Anatoly Sorokin and Igor Goryanin_

# Outline:

1. Summary
2. Introduction
3. Results
4. Parameters
5. Installation
6. Data Preparation
7. Manual download of files
8. How to cite us
9. References


## Summary
**What is it?** The functional and taxonomic analysis is the critical step in understanding the interspecies interaction within the biological community. At the moment these types of analysis are run separately which makes results difficult for interpretation. Here we present ASAR interactive tool for simultaneous analysis of the metagenomic data along three dimensions: taxonomy, function, metagenome.

**What is required from a user?** The user needs to have R and some R packages installed in the user’s computer. Additionally, terminal is needed to run the BASH script that will download, process and save data in the appropriate format to be read by the application. The annotation files from MG-RAST are used as an example data. In order to use data from MG-RAST only project ID and user webkey are needed from the user.

**Why should I use it?** Advantages of the tool are 1) Integrated functional and taxonomic analysis; 2) Comparative analysis of KEGG pathway enrichments; 3) KEGG Pathway Maps; 4) User-friendly interface.

## Introduction
This application implements two general analyses. First, it builds 3D dataset consisting of axes of taxonomy, functions and metagenome samples [1,2]. Second, through KEGG metabolic pathways analysis consists of a comparative analysis of pathways enrichment and visualization of pathways themselves [3].

Since taxonomic and functional annotations have many groups at several levels and metagenome samples are numerous, two main data manipulations are implemented. First, “Selection” involves selecting one or several separate groups in each axis at a certain level in functions and taxonomy. Second, “Aggregation” involves selecting a level lower than “Selection” level at which selected data should be aggregated to groups of this level by summing read counts. “Aggregation” of metagenomes is done by averaging metagenomes with same defined name.

There is the sample data to explore the app, and to use your own data, please, read the section “Data preparation”.

## Results

### _3D dataset (Function & Taxonomy & Metagenomes) Interactive Heatmaps_

Heatmaps contain dendrograms to the left and above from heatmap and column and row names below and to the right from the heatmap. Dendrograms are generated with “hclust” function in R and color key upleft to heatmaps represents the distribution of color in the heatmap. The value of the cell in the heatmap can be viewed by locating mouse cursor above that cell. The value of the cell is the log2 value of read count for that cell.

![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image1.png)
_Figure 1._ **Function vs Taxonomy (F/T)**

![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image002.gif)
_Figure 2._ **Function vs Metagenomes (F/M)**

![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image003.gif)
_Figure 3._ **Taxonomy vs Metagenomes (T/M)**

**KEGG Pathways Heatmap**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image004.gif)
_Figure 4._ **KEGG Pathways Heatmap**

In addition to characteristics explained for 3D dataset heatmaps, KEGG Pathways Heatmap has Standard Deviation cutoff explained in parameters section for SD cutoff.

**KEGG Pathway**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image005.gif)
_Figure 5._ **KEGG Pathway**

KEGG Pathway shows the metabolic pathway image for selected KEGG pathway and color enzymes with known expression in metagenomes. Each rectangle representing enzyme is horizontally divided into a number of selected metagenomes, where coloration of partitions corresponds to values of selected metagenomes(i.e. Partitions are colored from left to right as an order of selected metagenomes). Values of enrichment represent percent contribution of selected taxon from total read count of the whole metagenome to this KEGG Orthology(enzyme). The color key will represent a variation from zero to maximum value among KO’s in the pathway.

## Parameters
### _Metagenome Selection_

Function vs Taxonomy (F/T) heatmap requires selection of single metagenome, while all others allow selection of multiple metagenomes(see figure 6). Names of metagenomes shown for selection can be changed by selecting a specific column of metadata, rows of which will be used as metagenome names (see figure 7). In order to consider several metagenomes as one, these metagenomes should be given the same name in selected metadata column. In this case abundances of reads of metagenomes with same names are averaged in new metagenome under shared name before any abundances processing steps.

![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image006.gif)
_Figure 6._ **Unlike other heatmaps, F/T heatmap requires selection of single metagenome.**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image007.gif)
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image008.gif)

**Taxonomy and Function Selection**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image009.gif)
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image010.gif)

There are 8 choices for taxonomy level, where "root" means all domains and 7 levels from “domain” to “strain”. At all levels except "root" taxon selection will give a list of all available taxons at that level and multiple taxon selections are possible. At last, second taxon level is used to aggregate selected taxons at that level. Functions are selected with the same principle with the only difference that it has 5 levels instead of 8 (see figure 8).

**Pathway Selection for Building KEGG Pathway**

The list of KEGG Pathways available for current selection of metagenomes and taxons is displayed for selection. Selecting Pathway will send the request to build KEGG Pathway and may take up to several seconds depending on the number of genes in the pathway.

**Heatmap Height Selection**

Default height of heatmaps is 20 pixels per row and is adjustable through slider input parameter, which displays and sets the height of a single row in pixels.

**Image Download**

Every heatmap can be downloaded by typing a user-defined file name in “Enter file name” text input parameter and selecting image format (PNG or PDF) and subsequently pressing the Download button. KEGG map image is downloaded in the same way but without defining file name. Alternatively, map image can be saved by clicking the right mouse button (usual browser functionality), where user can define filename.

**Standard deviation cutoff for KEGG Orthology terms**

KEGG Pathways Heatmap and KEGG Pathway have the parameter called “SD cutoff for KO terms” which defines the value of Standard Deviation for individual KO’s among all selected metagenomes and is adjustable. This value cuts off all KO’s with SD less than that value from Heatmap.

**Parameters in Metadata Tab**

Selection of column of metadata displays all column of metadata and allows to use rows of these columns as names of metagenomes in metagenome selection parameter. Metadata values are editable and new columns with a name specified by user can be added to the metadata. There are three types of new column user can select, namely “integer”, “double” and “character”. Pushing “Save” button updates default RData file with current dataset. Pushing “Save a new column” button will update metadata but current Rdata should be saved and run in new session. When metagenome names from some column of metadata are used identically named metagenome data will be averaged and analyzed as single metagenome. This is how aggregation of metagenomes is done.

**Parameters in Settings Tab**

Upload function can be used to upload Rdata files generated previously and Save function can be used to save current state of loaded dataset. Settings tab has entries for changing default values of functional and taxonomic levels. Palette of colors used for coloring heatmaps can also be selected. Pressing “Save changes” button will save these default parameters for next sessions.

## Installation

R Packages from CRAN:

Package ‘shiny’ _version 1.0.3_
Package ‘ggplot2’ _version 2.2.1_
Package ‘gplots’ _version 3.0.1_
Package ‘data.table’ _version 1.10.4_
Package ‘plyr’ _version 1.8.4_
Package ‘stringr’ _version 1.2.0_
Package ‘shinythemes’ _version 1.1.1_
Package ‘matrixStats’ _version 0.52.2_
Package ‘png’ _version 0.1-7_
Package ‘devtools’ _version 1.13.2_
Package ‘rhandsontable’ _version 0.3.4.6_
Package ‘RColorBrewer’ _version 1.1-2_

Bioconductor:

Package ‘mmnet’ _version 1.15.0-1_
Package ‘pathview’ _version 1.14.0_
Package ‘biomformat’ _version 1.2.0_
Package ‘KEGGREST’ _version 1.14.1_
Package ‘limma’ _version 3.30.13_

GitHub:

Package ‘d3heatmap’ by _“Alanocallaghan/d3heatmap”_

**To run the app on your local machine:**

1.Download RStudio/R

2.Run these commands below in the console:

```markdown
install.packages(c("shiny","ggplot2","gplots","RColorBrewer","data.table","plyr","stringr","shinythemes","matrixStats","png","devtools","rhandsontable"))

##try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("mmnet")

biocLite("pathview")

biocLite("biomformat")

biocLite("KEGGREST")

biocLite("limma")

library(devtools)
install_github("Alanocallaghan/d3heatmap")

#Run this command in the console: 

shiny::runGitHub("ASAR", "Askarbek-orakov",subdir="R")

```

## Data Preparation

Preparation of data from MG-RAST only requires project ID and webkey to be given to BASH script which subsequently downloads and processes all required files and generates Rdata file that is directly used by the application.

## Manual Download of Files

_The list of required input files:_

1.Functional annotations file either by KEGG or SEED
2.Taxonomic annotations file either by KEGG or SEED
3.KEGG Orthology file
4.Biom file
5.Metadata

Our app uses MG-RAST annotations as an example. MG-RAST has both public and private projects which can be downloaded as it is described in its manual. There are two ways of downloading files from MG-RAST and prepare them for the input.

The first way is to download files directly from MG-RAST website, while second way is to download files through API or other command line tools, such as a terminal. In the case of the former way, you will have to rename files manually, while in the case of the latter way, given code will download files automatically and rename them automatically.




![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image011.gif)
_Figure 11._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image012.gif)
_Figure 12._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image013.gif)
_Figure 13._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image014.gif)
_Figure 14._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image015.gif)
_Figure 15._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image016.gif)
_Figure 16._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image017.gif)
_Figure 17._ **Function vs Metagenomes (F/M)**
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image018.gif)
_Figure 18._ **Function vs Metagenomes (F/M)**

```markdown
Syntax highlighted code block

# Nazgul
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/nazgul-sakenova01/test2/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
