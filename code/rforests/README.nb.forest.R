if ( F )
{
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
  source(paste(RHOME,"variation.filtering.R",sep="/"))
  source(paste(RHOME,"bsub.R",sep="/"))
  source("README.databases.R")
  library(ROC)
}
#do.classify <- F
#do.summarize <- T
#do.bsub <- T
if ( T )
{
  stub <- "xchip"
  fname <- paste("nb.forest.boot",
                 if(is.null(rdate)) gsub("-","",as.character(Sys.Date())) else rdate,
                 sep=".")
  cat( "Results will be saved to '", fname, "'\n", sep="" )
  
  SUMMARY <- matrix( NA, nrow(DBs), 5, dimnames=list(DBs[,"ID"],NULL) )
  BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
  names(BSUB.OUT) <- DBs[,"ID"]
  
  if (do.classify) {
    cat( "Results will be saved to '", fname, "'\n", sep="" )
  }
  if (do.summarize) {
    fname <- paste("nb.forest.boot",rdate,sep=".")
  }
  missing <- NULL
  #for ( i in match(c("rosetta.nejm","AllAmlMLL"),DBs[,"ID"]) )
  for ( i in DB.idx )
  {
    if (do.classify)
    {
      cat( "*****************************\n" )
      cat( "Processing", DBs[i,"ID"], "\n" )
      cat( "*****************************\n" )
      
      cmd <- paste(
                   "cat( \"NB-forest classification\\n\" )",
                   "library(ROC)",
                   "library(SparseM)",
                   "CLS <- read.cls( paste(stub,DBs[i,\"CLS\"],sep=\"/\") )",
                   "DAT <- read.data.gp( paste(stub,DBs[i,\"DB\"],sep=\"/\") )",
                   "if(nrow(DAT)>4000) DAT <- variation.filter( DAT, ngenes=4000 )",
                   "OUT <- nb.forest.estimate(dat=t(DAT@signal),cls=CLS,ntree=200,resample=\"boot\",exclude=\"none\",do.importance=T,nrnd=1000,verbose=T,score=\"ACC\",maxworse=3,max.nfeat=30)",
                   "save(OUT,file=paste(fname,DBs[i,\"ID\"],sep=\".\"))\n",
                   sep="\n")
      BSUB.OUT[[i]] <- bsub(cmd)
    }
    if (do.summarize)
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
      
      outI <- cls.summary(nb.forest.summary(OUT),CLS)[1:5] 
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
