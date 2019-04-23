if ( F ) {
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
if ( T )
{
  for ( resample in c("boot","nbf","nba") ) {
  #for ( resample in c("nbf","nba") ) {
  stub <- "xchip"
  fname <- paste("rpart.forest", resample,
                 if(is.null(rdate)) gsub("-","",as.character(Sys.Date())) else rdate,
                 sep=".")
  cat( "Results will be saved to '", fname, "'\n", sep="" )
  
  SUMMARY <- matrix( NA, nrow(DBs), 5, dimnames=list(DBs[,"ID"],NULL) )
  DETAILS <- BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
  names(DETAILS) <- names(BSUB.OUT) <- DBs[,"ID"]
  
  model <- list(estimate=rpart.forest.estimate,
                predict=rpart.forest.predict,
                predictors=rpart.forest.predictors)
  
  missing <- NULL
  #for ( i in match(c("dlbcl","duke","duke.er","lung","rosetta","rosetta.er"),DBs[,"ID"]) )
  for ( i in DB.idx )
  {
    if ( do.classify )
    {
      cat( "*****************************\n" )
      cat( "Processing", DBs[i,"ID"], "\n" )
      cat( "*****************************\n" )
      
      cmd <- paste(
                   paste("cat( \"Rpart (",resample,") ensemble classification\\n\" )",sep=""),
                   "library(ROC)",
                   "library(rpart)",
                   "library(survival)",
                   "CLS <- read.cls( paste(stub,DBs[i,\"CLS\"],sep=\"/\") )",
                   "DAT <- read.data.gp( paste(stub,DBs[i,\"DB\"],sep=\"/\") )",
                   "if(nrow(DAT)>3000) DAT <- variation.filter( DAT, ngenes=3000 )",
                   "DAT <- t(DAT@signal)",
                   "colnames(DAT) <- paste( \"G\", 1:ncol(DAT), sep=\"\")",
                   "nfolds <- if(nrow(DAT)>250) min(10,min(tabulate(CLS+1))-1) else nrow(DAT)",
                   #"prior <- rep(1/nlevels(as.factor(CLS)),nlevels(as.factor(CLS)))",
                   "prior <- tabulate(as.factor(CLS))/length(CLS)",
                          if ( resample=="boot" ) {
                     "OUT <- cv.classify(dat=DAT,cls=CLS,model=model,nfeats=ncol(DAT),ntree=50,nfolds=nfolds,resample=\"boot\",exclude=\"none\",oob=F,parms=list(prior=prior))"
                   }
                   else if ( resample=="nbf" ) {
                     "OUT <- cv.classify(dat=DAT,cls=CLS,model=model,nfeats=ncol(DAT),ntree=50,nfolds=nfolds,resample=\"none\",exclude=\"first\",oob=F,parms=list(prior=prior))"
                   }
                   else if ( resample=="nba" ) {
                     "OUT <- cv.classify(dat=DAT,cls=CLS,model=model,nfeats=ncol(DAT),ntree=50,nfolds=nfolds,resample=\"none\",exclude=\"all\",oob=F,parms=list(prior=prior))"
                   }
                   else {
                     stop( "resample type not recognized: ", resample )
                   },
                   "save(OUT,file=paste(fname,DBs[i,\"ID\"],sep=\".\"))\n",
                   sep="\n")
      BSUB.OUT[[i]] <- bsub(cmd,options="-M 2000" )
      
      cat( "\n" )
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
}
