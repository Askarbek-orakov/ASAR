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

1. Package ‘shiny’ _version 1.0.3_
2. Package ‘ggplot2’ _version 2.2.1_
3. Package ‘gplots’ _version 3.0.1_
4. Package ‘data.table’ _version 1.10.4_
5. Package ‘plyr’ _version 1.8.4_
6. Package ‘stringr’ _version 1.2.0_
7. Package ‘shinythemes’ _version 1.1.1_
8. Package ‘matrixStats’ _version 0.52.2_
9. Package ‘png’ _version 0.1-7_
10. Package ‘devtools’ _version 1.13.2_
11. Package ‘rhandsontable’ _version 0.3.4.6_
12. Package ‘RColorBrewer’ _version 1.1-2_

Bioconductor:

1. Package ‘mmnet’ _version 1.15.0-1_
2. Package ‘pathview’ _version 1.14.0_
3. Package ‘biomformat’ _version 1.2.0_
4. Package ‘KEGGREST’ _version 1.14.1_
5. Package ‘limma’ _version 3.30.13_

GitHub:

Package ‘d3heatmap’ by _“Alanocallaghan/d3heatmap”_


**To run the app on your local machine:**

1. Download RStudio/R

2. Run these commands below in the console:

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

1. Functional annotations file either by KEGG or SEED
2. Taxonomic annotations file either by KEGG or SEED
3. KEGG Orthology file
4. Biom file
5. Metadata

Our app uses MG-RAST annotations as an example. MG-RAST has both public and private projects which can be downloaded as it is described in its manual. There are two ways of downloading files from MG-RAST and prepare them for the input.

The first way is to download files directly from MG-RAST website, while second way is to download files through API or other command line tools, such as a terminal. In the case of the former way, you will have to rename files manually, while in the case of the latter way, given code will download files automatically and rename them automatically.


I  ### **Download files directly from MG-RAST website as follows.**

a. Select a project.

![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image011.gif)

b. Select a sample.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image012.gif)

c. Go to Download.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image013.gif)

d. Download functional annotation file either by SEED or KEGG by selecting “function” for **Annotation Type** and either “KEGG” or “SEED” for **Data Source** and rename them by adding “.fkegg” or “.fseed” respectively. 
_Examples:_ “mgm4714675.3.fkegg” and “mgm4714675.3.fseed”.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image014.gif)
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image015.gif)


e. Download taxonomic annotation file either by SEED or KEGG by selecting “organism” for **Annotation Type** and either “KEGG” or “SEED” for **Data Source** and rename them by adding “.kegg” or “.seed” respectively. 
_Examples:_ “mgm4714675.3.kegg” and “mgm4714675.3.seed”.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image016.gif)
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image017.gif)


f. Download KEGG Orthology file by selecting “ontology” for **Annotation Type** and “KO” for **Data Source** and rename the file by adding at the end “.ko”. 
_Example:_ “mgm4714675.3.ko”.

g. Biome file can be downloaded only from MG-RAST API/command line tools (see below)

h. Download metadata file by entering a project and pressing file icon as shown below and rename the file as “jobs.tsv”.

1.Press icon shown below.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image018.gif)

2. Select all.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image018.gif)

3. Press the icon shown below to download metadata.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image018.gif)

II  ### **Download files through API or other command line tools, such as a terminal.**

1. You can find how to download files through API by clicking this link.

2. An example of how to download files through terminal is shown below.

Webkey will be needed to download files from private projects. To get your webkey in MG-RAST, press “ show webkey” as indicated below.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image018.gif)

Open Terminal and run this chunk of code after modifying webkey and metagenome ID (as marked by red squares) to download your files.
![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image018.gif)

1. Functional annotation by SEED
```markdown
curl -H "auth: your_webkey_comes_here" -H 'Accept-Encoding: gzip,deflate' "http://api-pql.metagenomics.anl.gov/1/annotation/similarity/mgm4714675.3?source=SEED&type=function&identity=60&length=15" -o mgm4714679.3.fseed
```
2. Taxonomic annotation by SEED
```markdown
curl -H "auth: your_webkey_comes_here" -H 'Accept-Encoding: gzip,deflate' "http://api.metagenomics.anl.gov/1/annotation/similarity/mgm4714675.3?source=SEED&type=organism&identity=60&length=15" -o mgm4714675.3.seed
```
3. Functional annotation by KEGG
```markdown
curl -H "auth: your_webkey_comes_here" -H 'Accept-Encoding: gzip,deflate' "http://api-pql.metagenomics.anl.gov/1/annotation/similarity/mgm4714675.3?source=KEGG&type=function&identity=60&length=15" -o mgm4714679.3.fkegg
```
4. Taxonomic annotation by KEGG
```markdown
curl -H "auth: your_webkey_comes_here" -H 'Accept-Encoding: gzip,deflate' "http://api-pql.metagenomics.anl.gov/1/annotation/similarity/mgm4714675.3?source=KEGG&type=organism&identity=60&length=15" -o mgm4714679.3.kegg
```
5. KEGG Orthology file
```markdown
curl -H "auth: your_webkey_comes_here" -H 'Accept-Encoding: gzip,deflate' "http://api-pql.metagenomics.anl.gov/1/annotation/similarity/mgm4714675.3?source=KO&type=ontology&identity=60&length=15" -o mgm4714675.3.ko
```
6. Biome file:
```markdown
curl -H "auth: your_webkey_comes_here" -H 'Accept-Encoding: gzip,deflate' "http://api-pql.metagenomics.anl.gov/1/matrix/organism?id=mgm4714675.3&id=mgm4714661.3&id=mgm4714663.3&source=SEED&group_level=strain&result_type=abundance&hit_type=all&identity=60&length=15" -o mgm.biome
```
7. Metadata
Metadata is created by user.

**Creating .RData**

Put all five files into the same directory.
You may use this app by
A. Exploring the pre-loaded example data set. This is a pre-loaded Metagome Samples taken from swine waste example for exploring the app's features.

B. Upload your own data that is either
C. Uploading an .RData file containing your data that was previously downloaded from the app session.

## Cite us:

Askarbek N. Orakov, Nazgul K. Sakenova, Anatoly Sorokin and Igor Goryanin. 2017. “ASAR: ”.

## References:

[1] Keegan, K. P., Glass, E. M., & Meyer, F. (2016). MG-RAST, a metagenomics service for analysis of microbial community structure and function. Microbial Environmental Genomics (MEG), 207-233.

[2] Overbeek, R., Begley, T., Butler, R. M., Choudhuri, J. V., Chuang, H. Y., Cohoon, M., ... & Fonstein, M. (2005). The subsystems approach to genome annotation and its use in the project to annotate 1000 genomes. Nucleic acids research, 33(17), 5691-5702.

[3] Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., & Tanabe, M. (2016). KEGG as a reference resource for gene and protein annotation. Nucleic acids research, 44(D1), D457-D462.


For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/nazgul-sakenova01/test2/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
