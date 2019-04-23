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
  nfeats <- c( 2,5,10,20,30,50,100,500,1000,0 )
  fname <- paste("naive.cv2",
                 if(is.null(rdate)) gsub("-","",as.character(Sys.Date())) else rdate,
                 sep=".")
  
  SUMMARY <- matrix( NA, nrow(DBs), 5, dimnames=list(DBs[,"ID"],NULL) )
  DETAILS <- BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
  names(DETAILS) <- names(BSUB.OUT) <- DBs[,"ID"]
  
  if (do.classify) {
    cat( "Results will be saved to '", fname, "'\n", sep="" )
  }
  if (do.summarize) {
    fname <- paste("naive.cv2",rdate,sep=".")
  }
  missing <- NULL
  #for ( i in match(c("rosetta.nejm","AllAmlMLL"),DBs[,"ID"]) )
  for ( i in DB.idx )
  {
    if ( do.classify ) {
      cat( "*********************************************\n" )
      cat( "Processing ", DBs[i,"ID"], " (", DBs[i,"DB"], ")\n", sep="" )
      cat( "*********************************************\n" )

      cmd <- paste(
                   "cat( \"Naive-cv2 classification\\n\" )",
                   "library(ROC)",
                   "DAT <- read.data.gp( paste(stub,DBs[i,\"DB\"],sep=\"/\") )@signal",
                   "CLS <- read.cls( paste(stub,DBs[i,\"CLS\"],sep=\"/\") )",
                   "nlevs <- length(levels(CLS))",
                   "nfeatsI <- nfeats*nlevs; nfeatsI[nfeatsI>nrow(DAT)] <- 0; nfeatsI <- unique(nfeatsI)",
                   "cat( \"features tested:\", paste(nfeatsI,sep=\",\"), \"\\n\" )",
                   "nfolds2 <- min(10,min(tabulate(CLS+1))-1)",
                   "OUT <- naive.cv2(dat=t(DAT),cls=CLS,nfeats=nfeatsI,nfolds1=ncol(DAT),nfolds2=nfolds2,details=T,verbose=T)\nsave(OUT,file=paste(fname,DBs[i,\"ID\"],sep=\".\"))\n",
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
      
      outI <- cls.summary(OUT$probs,CLS)[1:5]
      SUMMARY[i,] <- outI
      if (i==1) colnames(SUMMARY) <- names(outI)[1:ncol(SUMMARY)]
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
