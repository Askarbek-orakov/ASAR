Askarbek Orakov, Nazgul Sakenova, Anatoly Sorokin and Igor Goryanin

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

![](https://github.com/Askarbek-orakov/ASAR/blob/master/docs/media/image001.gif)


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
