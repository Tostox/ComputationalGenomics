RHOME <- "~/dvlp/R"
library(MASS)
source(paste(RHOME,"heatmap.R",sep="/"))
source(paste(RHOME,"mysource.R",sep="/"))
source(paste(RHOME,"read.res.new.R",sep="/"))

sim.quadnorm <- function( M, N, Mu0, Mu1, Sigma, seed=NULL )
{
  if(!is.null(seed))
    set.seed(seed)

  DAT <- matrix(0,M,N)
  base <- offset <- 0
  for ( i in 1:round(M/4/4) ) # class 1 markers
  {
    base <- (i-1)*4 + 1
    DAT[base:(base+3),]<- cbind(t(mvrnorm(N/2,Mu0,Sigma)),
                                t(mvrnorm(N/2,Mu1,Sigma)))
  }
  #cat( "end:", base+3, "\n" )
  base <- base+3
  for ( i in 1:round(M/4/4) ) # class 0 markers
  {
    offset <- (i-1)*4 + 1
    DAT[(base+offset):(base+offset+3),] <- cbind(t(mvrnorm(N/2,Mu1,Sigma)),
                                                 t(mvrnorm(N/2,Mu0,Sigma)))
  }
  #cat( "end:", base+offset+3, "\n" )
  base <- base+offset+3
  for ( i in 1:((M-round(M/4/4)*8)/4) ) # random genes
  {
    offset <- (i-1)*4 + 1
    DAT[(base+offset):(base+offset+3),] <- t(mvrnorm(N,Mu0,Sigma))
  }
  #cat( "end:", base+offset+3, "\n" )
  colnames(DAT) <- paste("S",1:N,sep="")
  rownames(DAT) <- paste("G",1:M,sep="")
  DAT
}
M <- 1000
N <- 100
Mu0 <- c(0,0,0,0)
Mu1 <- c(0,0,4,4)
Sigma1 <- matrix( c(2,1,1,1,1,2,2,2,1,2,5,5,1,2,5,8), 4, 4)

DAT<- sim.quadnorm(M=M,N=N,Mu0=Mu0,Mu1=Mu1,Sigma=Sigma1,seed=12345)
DAT4 <- new("gctdata", signal=DAT, description=rownames(DAT))
CLS <- rep(0:1,each=N/2); levels(CLS) <- c("0","1")

#my.heatmap(DAT@signal,scale="r",ColSideColors=rep(c("green","orange"),each=N/2),RowSideColors=rep(c("orange","green","light blue"),times=c(M/4,M/4,M/2)))

if ( F )
{
  RF2 <- randomForest( x=t(DAT2@signal), y=as.factor(CLS), ntree=1000, importance=T )
  ord <- order(RF2$importance[,3],decreasing=T)
  zeros <- c(seq(3,500,4),seq(4,500,4))
  ones <- c(seq(1,497,4),seq(2,498,4))
  impOrd <- order(-RF2$importance[,3])
  importance <- rep(0,1000)
  importance[impOrd[1:500]] <- 1
  markers <- rep(1:0,each=500)
  ftable(importance,markers)

  markers <- rep(0,1000); markers[zeros] <- 1
  ftable(importance,markers)
  
  univar2 <- diffanal(DAT2,CLS,nperm=100,score=t.score)
  difOrd <- order(univar2$diff[,"asymp.p"])
  different <- rep(0,1000)
  different[difOrd[1:500]] <- 1
  markers <- rep(1:0,each=500)
  ftable(different,markers)
