## various pacakges (more than you probably need)

pkgs <- c("mclust",
          "getopt",
          "optparse",
          "randomForest",
          "randomForestSRC",
          "pamr",
          "e1071",
          "combinat",
          "heatmap.plus",
          "rgl",
          "ggplot2",
          "ggdendro",
          "gridExtra",
          "NanoStringNorm",
          "gplots",
          "rjags",
          "PCIT",
          "rjson",
          "xlsx",
          "dynamicTreeCut",
          "cba",
          "devtools",
          "roxygen2",
          "rmarkdown",
          "XLConnect",
          "data.table",
          "dendextend",
          "VennDiagram",
          "caret",
          "pROC",
          "klaR")
install.packages(pkgs,repos="http://cran.r-project.org")

## packages not available through the install.packages command
require(devtools)
install_github('hadley/staticdocs')
install_github("andrie/ggdendro")
install_github("hadley/ggplot2")

## install bioconductor and needed packages

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("Biobase","biomaRt","ASSIGN","pathifier","oligo","oligoClasses"))
biocLite(c("limma","frma","GSEAlm","ConsensusClusterPlus"))
biocLite(c("edgeR","DESeq2"))

## installation of CBMRtools (our own scripts)

#installing current tag
require(devtools)
PAT <- "04fe676593e46b6bda5a5d09431156e8a500349a"
install_github("montilab/CBMRtools/CBMRtools",ref="v1.0.3", auth_token = PAT)   
require(CBMRtools)

