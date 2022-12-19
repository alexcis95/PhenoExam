# PhenoExam

## How to install PhenoExam package

~~~
# If you don't have devtools package
install.packages("devtools")
library(devtools)

# Install PhenoExam package 
install_github("alexcis95/PhenoExam")
library(PhenoExam)
~~~

You need the following bioconductor packages:  

clusterProfiler  
AnnotationDbi  
org.Hs.eg.db  

~~~
# How to install Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
~~~

You need the following R packages:  

readr  
data.table  
ggplot2  
stats  
plotly  
ggpubr  
dplyr  
viridis 
parallel  
purrr  
DT  
Hmisc  
pheatmap  


## Uses of PhenoExam


**PhenoExam** is an R package that performs (1) phenotype enrichment analysis on a gene set, (2) measures statistically significant phenotype similarities between gene sets and (3) detects significant differential phenotypes for them. Phenotypic Similarity between two groups of genes is performed by assessing the statistical significance of the Phenotypic Overlap Ratio (POR) between those (i.e. the number of common relevant phenotypes between the gene sets). PhenoExam uses the HPO, MGD, and CRISPRbrain databases for phenotypes. PhenoExam uses UNIPROT, CTD, ORPHANET, CLINGE, GENOMICS ENGLAND, CGI and PSYGENET to obtain diseases terms through DisGeNET database.   


# Available tutorials

[Phenotype analysis tutorial](https://raw.githack.com/alexcis95/PhenoExamWebTutorials/main/tutorial.html)


# PhenoExamWeb Shiny app  

Instead using the shiny app [PhenoExamWeb shiny app](https://alejandrocisterna.shinyapps.io/phenoexamweb) it might be a better option to consider using the R package (https://github.com/alexcis95/PhenoExam) or the shiny app locally (https://github.com/alexcis95/PhenoExam/blob/master/PhenoExamWeb.zip) for heavy analyses. We recommend not to use the shiny app remotely if users want to perform a “comparator phenotype analysis” and the random test number (>40) because the shiny app server does not have enough RAM memory to run the analysis. 

If you want to use the shiny app in a local execution enviroment please download the [PhenoExamWeb.zip file](https://github.com/alexcis95/PhenoExam/blob/master/PhenoExamWeb.zip) and run the .rmd file in your machine. 


# About  


This work is the result of a pre-doctoral training contract for research staff financed by the CARM Ministry of Employment, Research and Universities, through the Seneca Foundation - Science and Technology Agency of the Region of Murcia. This package has developed by **Alejandro Cisterna García** as part of his PhD mentored by professor **Juan Antonio Botía Blaya**.  


Authors:
Alejandro Cisterna, Aurora González-Vidal, Daniel Ruiz, Jordi Ortiz, Alicia Gómez-Pascual, Zhongbo Chen, Mike Nalls, Faraz Faghri, John Hardy, Irene Díez, Paolo Maietta, Sara Álvarez, Mina Ryten, Juan A. Botía 
