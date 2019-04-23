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
do.bsub <- F

stub <- "xchip"
nfeats <- c( 2,5,10,20,30,50,100,500,1000,0 )
fname <- paste("naive.cv2",gsub("-","",as.character(Sys.Date())),sep=".")
cat( "Results will be saved to '", fname, "'\n", sep="" )

SUMMARY <- matrix( NA, nrow(DBs), 5, dimnames=list(DBs[,"ID"],NULL) )
DETAILS <- BSUB.OUT <- lapply(1:nrow(DBs), function(z) NULL )
names(DETAILS) <- names(BSUB.OUT) <- DBs[,"ID"]

for ( i in 1:nrow(DBs) )
{
  cat( "*****************************\n" )
  cat( "Processing", DBs[i,"ID"], "\n" )
  cat( "*****************************\n" )

  DAT <- read.data.gp( paste(stub,DBs[i,"DB"],sep="/") )@signal
  CLS <- read.cls( paste(stub,DBs[i,"CLS"],sep="/") )
  nlevs <- length(levels(CLS))
  nfeatsI <- nfeats*nlevs; nfeats <- nfeats[nfeats<nrow(DAT)]
  cat( "features tested:", paste(nfeatsI,sep=","), "\n" )

  if ( do.bsub )
  {
    cmd <- "OUT <- naive.cv2(dat=t(DAT),cls=CLS,nfeats=nfeatsI,nfolds1=ncol(DAT),nfolds2=10,details=T,verbose=T)\nsave(OUT,file=paste(fname,DBs[i,\"ID\"],sep=\".\"))\n"
    BSUB.OUT[[i]] <- bsub(cmd)
  }
  else
  {
    OUT <- naive.cv2(dat=t(DAT),cls=CLS,nfeats=nfeatsI,nfolds1=ncol(DAT),nfolds2=10,
                   details=T,verbose=T)
    DETAILS[[i]] <- OUT
    
    OUT <- cls.summary(OUT$probs,CLS)[1:5]
    SUMMARY[i,] <- OUT
    if (i==1) colnames(SUMMARY) <- names(OUT)[1:ncol(SUM)]
    save(SUMMARY,DETAILS,file=fname)
  }
  cat( "\n" )
}
if ( do.summary && do.bsub )
{
  for ( i in 1:nrow(DBs) )
  {
    cat( "*****************************\n" )
    cat( "Processing", DBs[i,"ID"], "\n" )
    cat( "*****************************\n" )
    
    load( file=paste(fname,DBs[i,"ID"],sep=".")
  }
}
