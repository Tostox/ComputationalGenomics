RHOME <- "~/dvlp/R"
source(paste(RHOME,"mysource.R",sep="/"))
source(paste(RHOME,"misc.math.R",sep="/"))
source(paste(RHOME,"read.res.new.R",sep="/"))  
source(paste(RHOME,"ftable.R",sep="/"))  
source(paste(RHOME,"gene.cluster.R",sep="/"))
source(paste(RHOME,"variation.filtering.R",sep="/"))
source(paste(RHOME,"bsub.R",sep="/"))
source("README.databases.R")
library(ROC)
stub <- "xchip"

BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
names(BSUB.OUT) <- DBs[,"ID"]
do.diffanal <- F
do.importance <- F
do.compare <- T
do.bsub <- F

if (do.diffanal)
{
  #for ( i in 1:nrow(DBs) )
  for ( i in match(c("dlbcl"),DBs[,"ID"]) )
  {
    command <- paste("bsub -o bsub.out.diffanal.", DBs[i,"ID"], " diffanal ",
                     stub, "/", DBs[i,"DB"], " ",
                     stub, "/", DBs[i,"CLS"],
                     " nperm=25000 output=", stub, "/", DBs[i,"ID"],".diffanal.xls ",
                     DBs[i,"PARAM"], " -exhaustive -verbose", sep="")
    cat( i, ") ", command, "\n", sep="" )
    system( command )
  }
}
if ( do.importance )
{
  fname <- "randomForest.importance"
  for ( mtry in c(2,3) )
  {
    SUMMARY <- matrix( NA, nrow(DBs), 2, dimnames=list(DBs[,"ID"],c("ERR","RERR")))
    fnameStub <- paste(fname,".mtry", mtry, sep="")
      
    for ( i in 1:nrow(DBs) )
    #for ( i in match(c("rosetta.nejm"),DBs[,"ID"]) )
    {
      cat( "*****************************\n" )
      cat( "Processing", DBs[i,"ID"],   "\n" )
      cat( "*****************************\n" )
      
      fnameO <- paste(fnameStub, DBs[i,"ID"],"xls",sep=".")
      fnameR <- paste(fnameStub, DBs[i,"ID"],"R.out",sep=".")
      
      cmd <- paste("library(randomForest)",
                   "CLS <- read.cls( paste(stub,DBs[i,\"CLS\"],sep=\"/\") )",
                   "DAT <- read.data.gp( paste(stub,DBs[i,\"DB\"],sep=\"/\") )",
#                  "if(nrow(DAT)>3000) DAT <- variation.filter( DAT, ngenes=3000 )",
                   "DAT <- t(DAT@signal)",
                   "colnames(DAT) <- paste( \"G\", 1:ncol(DAT), sep=\"\")",
                   "OUT <- randomForest(x=as.data.frame(DAT),y=as.factor(CLS),ntree=1000,mtry=mtry*floor(sqrt(ncol(DAT))),importance=T,keep.forest=F)",
                   "my.write.matrix(signif(OUT$confusion,4),sep=\"\\t\",row.names=T,file=fnameO)",
                   "cat(\"\\n\",append=T,file=fnameO)",
                   "my.write.matrix(signif(OUT$importance,3),sep=\"\\t\",row.names=T,append=T,file=fnameO)",
                   "save(OUT,file=fnameR)\n",
                   sep="\n")
      
      if (do.bsub) {
        BSUB.OUT[[i]] <- bsub(cmd,options="-M 2000 -J RFimp")
        next
      }
      load(file=fnameR)
      confusion <- OUT$confusion[,-ncol(OUT$confusion)]
      SUMMARY[i,] <- c(ERR=1 - sum(diag(confusion))/sum(confusion),
                       RERR=mean(1-diag(confusion)/apply(confusion,1,sum)))
    }
    if (!do.bsub) {
      cat( "Writing summary to ", fnameStub, ".sum.xls .. ", sep="" )
      my.write.matrix(signif(SUMMARY,4),sep="\t",row.names=T,file=paste(fnameStub,".sum.xls",sep=""))
      cat( "done.\n" )
    }
  }  
}
if( do.compare )
  pdf(paste("randomForest.importance.mtry",2,".png",sep=""))
if ( do.compare ) for ( i in 1:nrow(DBs))
{
  difName <- paste(stub, "/", DBs[i,"ID"],".diffanal.xls", sep="")
  impName <- paste("randomForest.importance.mtry",2,".",DBs[i,"ID"],".R.out",sep="")
  datName <- paste(stub,DBs[i,"DB"],sep="/")
  if (file.access(difName)[1]!=0 ||
      file.access(impName)[1]!=0 ||
      file.access(datName)[1]!=0) next
  
  cat( "*****************************\n" )
  cat( "Processing", DBs[i,"ID"],   "\n" )
  cat( "*****************************\n" )
  
  load(file=impName)
  if ( length(levels(OUT$y))>2 )
    next
  DIF <- read.delim( difName, header=T, row.names=1, skip=1 )
  DAT <- read.data.gp(datName)
  importance <- OUT$importance
  rownames(importance) <- genenames(DAT)
  idxI <- match(rownames(DIF),rownames(importance))
  if (any(is.na(idxI))) stop( "importance and dif don't match")
  importance <- importance[idxI,]
  N <- 200
  topNimp <- order(importance[,3],decreasing=T)[1:N]
  plotName <- paste("randomForest.importance.mtry",2,".",DBs[i,"ID"],".png",sep="")
  cat( "plotting ", plotName, " .. ", sep="" )
  plot(DIF[,"p2"],importance[,3],pch=".",xlab="diffanal p-value",ylab="importance (MeanDecreaseAccuracy)",
       main=DBs[i,"ID"])
  points(DIF[topNimp,"p2"],importance[topNimp,3],pch="*",col="red")
  cat( "done.\n" )
}
dev.off()
