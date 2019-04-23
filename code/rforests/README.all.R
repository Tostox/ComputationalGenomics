RHOME <- "~/dvlp/R"
source(paste(RHOME,"mysource.R",sep="/"))
source(paste(RHOME,"sample.R",sep="/"))
source(paste(RHOME,"naive.R",sep="/"))
source(paste(RHOME,"nb.forest.R",sep="/"))
source(paste(RHOME,"rpart.forest.R",sep="/"))
source(paste(RHOME,"misc.math.R",sep="/"))
source(paste(RHOME,"read.res.new.R",sep="/"))  
source(paste(RHOME,"xval.R",sep="/"))  
source(paste(RHOME,"cv.classifier.R",sep="/"))  
source(paste(RHOME,"fs.ova.R",sep="/"))  
source(paste(RHOME,"ftable.R",sep="/"))  
source(paste(RHOME,"gene.cluster.R",sep="/"))
source(paste(RHOME,"variation.filtering.R",sep="/"))
source(paste(RHOME,"bsub.R",sep="/"))
source("README.databases.R")
library(ROC)

do.classify <- F
do.summarize <- T
do.bsub <- F
rdate <- "20061031"
rdate <- "20060801"

do.classify <- T
do.summarize <- F
do.bsub <- T

DB.idx <- 1:nrow(DBs)
DB.idx <- match( c("mlbcl","egfr"), DBs[,"ID"] )
DB.idx <- match( "egfr", DBs[,"ID"] )
DB.idx <- match( c("duke.er"), DBs[,"ID"] )
DB.idx <- match( c("rosetta", "AllAmlMLL", "dmap13"), DBs[,"ID"] )

source( "README.naive.cv2.R" )
source( "README.stepwise.nb.R" )
source( "README.nb.forest.R" )
source( "README.rpart.R" )
source( "README.rpart.forest.R")
source( "README.randomForest.R" )

CNAMES <- c("naive.cv2",
            "stepwise.nb",
            "nb.forest.boot",
            "rpart",
            "rpart.forest.boot",
            "rpart.forest.nbf",
            "rpart.forest.nba",
            "randomForest")

SUMMARY <- matrix(NA, length(CNAMES)*2, nrow(DBs),
                  dimnames=list(as.vector(sapply(CNAMES,paste,c("ERR","RERR"),sep=".")),
                                DBs[,"ID"]))

for ( stub in CNAMES )
{
  cat( "Processing", stub, ".." )
  fname <- paste( stub, rdate, "xls", sep="." )
  sumry <- read.table( fname, sep="\t", header=T )
  idx <- grep( paste(stub,"[R]*ERR",sep="."), rownames(SUMMARY), fixed=F, perl=T )
  cat( "(idx: ", paste(idx,collapse=","), ") ", sep="" )
  SUMMARY[idx,] <- t(sumry[,c("ERR","RERR")])
  cat( "done.\n" )
}
my.write.matrix(round(SUMMARY*100,2),file=paste("summary.",rdate, ".xls",sep="."),
                row.names=T,names="DB")

outcome <- match(c("duke","dlbcl","prostate","rosetta","rosetta.nejm"),DBs[,"ID"])
plot(SUMMARY["nb.forest.boot.ERR",],SUMMARY["stepwise.nb.ERR",],pch=19,xlab="NB-forest",ylab="Stepwise-NB",main="Error rate comparison")
points(SUMMARY["nb.forest.boot.ERR",outcome],SUMMARY["stepwise.nb.ERR",outcome],pch=19,col="green")
abline(0,1)
