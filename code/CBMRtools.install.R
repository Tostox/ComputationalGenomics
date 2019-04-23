## Set of simple instructions to install our own CBMRtools package, a
## set of scripts (under development) for the analysis of (mainly)
## gene expression data.

## if not installed already, you might need to install these bioconductor packages first

source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("oligo")
biocLite("oligoClasses")

## then, do the following

library(devtools)
PAT <- "04fe676593e46b6bda5a5d09431156e8a500349a"
install_github("montilab/CBMRtools/CBMRtools",ref="v1.0.1", auth_token = PAT)

## the commands above are only needed once. Thereafter, you only need
## to include the following command

require(CBMRtools)
