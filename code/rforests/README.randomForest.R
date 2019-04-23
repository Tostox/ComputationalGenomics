if ( F )
{
  RHOME <- "~/dvlp/R"
  source(paste(RHOME,"mysource.R",sep="/"))
  source(paste(RHOME,"naive.R",sep="/"))
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
}
#do.bsub <- T
#do.classify <- F
#do.summarize <- T
if ( T )
{
  stub <- "xchip"
  fname <- paste("randomForest",
                 if(is.null(rdate)) gsub("-","",as.character(Sys.Date())) else rdate,
                 sep=".")
  cat( "Results will be saved to '", fname, "'\n", sep="" )
  
  SUMMARY <- matrix( NA, nrow(DBs), 5, dimnames=list(DBs[,"ID"],c("ERR","RERR","ROC","LS","NERR")) )
  BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
  names(BSUB.OUT) <- DBs[,"ID"]
  
  model <- list(
   estimate=function(dat,cls,...) randomForest(x=as.data.frame(dat),y=as.factor(cls),...),
   predict=function(dat,model) predict(object=model,newdata=as.data.frame(dat),type="prob"),
   predictors=function(obj) rownames(obj$importance[obj$importance[,"MeanDecreaseAccuracy"]>0,])
  )
  if (do.classify) {
    cat( "Results will be saved to '", fname, "'\n", sep="" )
  }
  missing <- NULL
  #for ( i in match(c("duke.er","dlbcl","prostate","rosetta","rosetta.er"),DBs[,"ID"]) )
  for ( i in DB.idx )
  {
    fnameI <- paste(fname,DBs[i,"ID"],sep=".")
    
    if ( do.classify )
    {
      cat( "*****************************\n" )
      cat( "Processing", DBs[i,"ID"], "\n" )
      cat( "*****************************\n" )
      
      cmd <- paste("cat( \"RandomForest classification\\n\" )",
                   "library(ROC)",
                   "library(randomForest)",
                   "CLS <- read.cls( paste(stub,DBs[i,\"CLS\"],sep=\"/\") )",
                   "DAT <- read.data.gp( paste(stub,DBs[i,\"DB\"],sep=\"/\") )",
                   "if(nrow(DAT)>4000) DAT <- variation.filter( DAT, ngenes=4000 )",
                   "DAT <- t(DAT@signal)",
                   "colnames(DAT) <- paste( \"G\", 1:ncol(DAT), sep=\"\")",
                   "OUT <- randomForest(x=DAT,y=as.factor(CLS),ntree=2000,importance=T,strata=as.factor(CLS))",
#                   "OUT <- cv.classify(dat=DAT,cls=CLS,model=model,nfeats=ncol(DAT),ntree=1000,importance=T)",
                   "save(OUT,file=fnameI)\n",
                   sep="\n")
      BSUB.OUT[[i]] <- bsub(cmd,options="-M 2000")
    }
    if ( do.summarize )
    {
      if (file.access(fnameI)[1]!=0) {
        missing <- c(missing, fnameI)
        next
      }      
      cat( "*********************************************\n" )
      cat( "Processing ", fnameI, "\n", sep="" )
      cat( "*********************************************\n" )
      
      load(file=fnameI)
      CLS <- read.cls( paste(stub,DBs[i,"CLS"],sep="/") )

      tmp <- OUT$confusion[,-ncol(OUT$confusion)]
      SUMMARY[i,"ERR"] <- (sum(tmp)-sum(diag(tmp)))/sum(tmp)
      SUMMARY[i,"RERR"] <- mean(OUT$confusion[,"class.error"])
      SUMMARY[i,"NERR"] <- sum(tmp)-sum(diag(tmp))
      #outI <- OUT$summary[nrow(OUT$summary),2:6]
      #SUMMARY[i,] <- outI
      #if (i==1) colnames(SUMMARY) <- names(outI)
    }
  }
  if (do.summarize) {
    fnameO <- paste(fname,".xls",sep="")
    cat( "Writing summary to '", fnameO, "' .. ", sep="" )
    my.write.matrix(round(SUMMARY,4),row.names=T,sep="\t",file=fnameO)
    cat( "done.\n" )
    if ( length(missing)>0 ) {
      cat( "files skipped:\n" )
      print( missing )
      cat( "\n" )
    }
  }
}
