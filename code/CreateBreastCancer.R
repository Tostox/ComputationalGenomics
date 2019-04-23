setwd("~/Research/Presentations/2015_10_09_BS830/")

## Installing required libraries
source("http://www.bioconductor.org/biocLite.R")
## in case you already have an installed version, but not the newest one:
## biocLite("BiocUpgrade")
biocLite("ArrayExpress","frma","limma","GSEAlm","affy")

## Getting the Richardson Breast Cancer dataset

## The dataset is available on Array Express under:
## http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-3744/

## However, there is a direct way to access the files from within R:
library(ArrayExpress)
AEset <- ArrayExpress("E-GEOD-3744",path='data',save=T)
colnames(pData(AEset)) <- gsub("Factor.Value..DiseaseState.","DiseaseState",colnames(pData(AEset)))
save(AEset,file='data/AEset.RData')
phenoTypeData<-pData(AEset)
write.table(phenoTypeData,file='data/phenotype.txt',sep='\t')
##load('rawdata/AEset.RData')

## fRMA Normalization
##
require(Biobase)
library(frma)
library(affy)
eSet <- frma(AEset, summarize = "robust_weighted_average")
saveRDS(eSet,file='data/normalizedBreastDB.RDS')
eSet <- readRDS('data/normalizedBreastDB.RDS')

## convert to gene symbols (does not work, using the one in 2013_09_06_BS830/code/analysis.Rmd)
##
require(biomaRt)
CBMDEV <- Sys.getenv('CBMDEV')
if (CBMDEV=="") stop( "Use 'setenv CBMDEV ..' to set CBMrepository's base directory" )
source( paste(CBMDEV, "R/misc.R", sep="/") )
source( paste(CBMDEV, "R/probesetAnnotation.R", sep="/") )
## check the filters by listFilters(mart)
eSet1 <- probesetAnnotation(eSet,filters="affy_hg_u133_plus_2")
saveRDS(eSet1,file='data/renamedBreastDB.RDS')
#eSet1 <- readRDS(file='data/renamedBreastDB.RDS') # from 2013_09_06_

eSet2 <- eSet1[,pData(eSet1)$DiseaseState %in%
               c("non-basal-like breast cancer","sporadic basal-like breast cancer")]
eSet2 <- eSet2[,order(pData(eSet2)$DiseaseState)]
    
## write GCT file
require(CBMRtools)

GCT <- new("gctdata",signal=exprs(eSet2),description=featureNames(eSet2))
CLS <- match(pData(eSet2)$DiseaseState,unique(pData(eSet2)$DiseaseState))-1
levels(CLS) <- unique(pData(eSet2)$DiseaseState)
write.gct(GCT,file='data/renamedBreastDB.gct')
write.cls(CLS,file='data/renamedBreastDB.cls')
