if (F) {
  RHOME <- "~/dvlp/R"
  source(paste(RHOME,"mysource.R",sep="/"))
  source(paste(RHOME,"naive.R",sep="/"))
  source(paste(RHOME,"nb.forest.R",sep="/"))
  source(paste(RHOME,"misc.math.R",sep="/"))
  source(paste(RHOME,"read.res.new.R",sep="/"))  
  source(paste(RHOME,"xval.R",sep="/"))  
  source(paste(RHOME,"cv.classifier.R",sep="/"))  
  source(paste(RHOME,"fs.ova.R",sep="/"))  
  source(paste(RHOME,"ftable.R",sep="/"))  
  source(paste(RHOME,"gene.cluster.R",sep="/"))
  source(paste(RHOME,"bsub.R",sep="/"))
  source("README.databases.R")
}
#do.classify <- T
#do.summarize <- F
#do.bsub <- T

if ( T )
{
  stub <- "xchip"
  nfeats <- c( 200 )
  fname <- paste("stepwise.nb",
                 if(is.null(rdate)) gsub("-","",as.character(Sys.Date())) else rdate,
                 sep=".")
  
  SUMMARY <- matrix( NA, nrow(DBs), 5, dimnames=list(DBs[,"ID"],NULL) )
  BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
  names(BSUB.OUT) <- DBs[,"ID"]
  
  model <- list(estimate=stepwise.nb,
                predict=naive.predict,
                predictors=names)
  
  if (do.classify) {
    cat( "Results will be saved to '", fname, "'\n", sep="" )
  }
  if (do.summarize) {
    fname <- paste("stepwise.nb",rdate,sep=".")
  }
  missing <- NULL
  #for ( i in match(c("leuk","rosetta.nejm","AllAmlMLL"),DBs[,"ID"]) )
  for ( i in DB.idx )
  {
    if (do.classify) {
      cat( "*****************************\n" )
      cat( "Processing", DBs[i,"ID"], "\n" )
      cat( "*****************************\n" )
      
      cmd <- paste(
                   "cat( \"Stepwise NB classification\\n\" )",
                   "library(SparseM)",
                   "library(ROC)",
                   "DAT <- read.data.gp( paste(stub,DBs[i,\"DB\"],sep=\"/\") )@signal",
                   "CLS <- read.cls( paste(stub,DBs[i,\"CLS\"],sep=\"/\") )",
                   "nlevs <- length(levels(CLS))",
                   "nfeatsI <- c( min(3000,nrow(DAT)) )",
                   "OUT <- cv.classify(dat=t(DAT),cls=CLS,model=model,nfeats=nfeatsI,max.nfeat=40,score=\"ACC\",do.bic=F,augmented=F,maxworse=5,nfolds=ncol(DAT))",
                   "save(OUT,file=paste(fname,DBs[i,\"ID\"],sep=\".\"))\n",
                   sep="\n")
      BSUB.OUT[[i]] <- bsub(cmd)
    }
    if ( do.summarize )
    {
      fnameI <- paste(fname,DBs[i,"ID"],sep=".")
      if (file.access(fnameI)[1]!=0) {
        missing <- c(missing, fnameI)
        next
      }
      cat( "*********************************************\n" )
      cat( "Processing ", fnameI, "\n", sep="" )
      cat( "*********************************************\n" )
      
      load(file=fnameI)
      CLS <- read.cls( paste(stub,DBs[i,"CLS"],sep="/") )
      
      outI <- OUT$summary[nrow(OUT$summary),2:6]
      SUMMARY[i,] <- outI
      if (i==1) colnames(SUMMARY) <- names(outI)
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
