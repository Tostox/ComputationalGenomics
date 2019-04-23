RHOME <- "~/dvlp/R"
source(paste(RHOME,"mysource.R",sep="/"))
source(paste(RHOME,"read.res.new.R",sep="/"))
source(paste(RHOME,"misc.math.R",sep="/"))
source(paste(RHOME,"xval.R",sep="/"))
source(paste(RHOME,"rpart.forest.R",sep="/"))
library(rpart)

surv <- read.table("/xchip/projects/confounders/datasets/rosetta_breast_out.annotation.txt",sep="\t",header=T)
surv <- Surv(surv[,"surtime"],surv[,"status"])
boxplot(surv[,1]~surv[,2],xlab="Status",ylab="Survival time (months)",col=c("dark green","dark red"))

dat <- read.res( "xchip/rosetta_breast_out.10Kmad.res" )

MAD <- apply(dat@signal,1,mad)
datN <- subset.resdata( dat, genenames=genenames(dat)[order(-MAD)[1:1000]] )
                                 
hatN <- rpart.forest.estimate( t(datN@signal), cls=surv, ntree=10 )
feats <- hatN$features[apply(hatN$features,1,sum)>0,]
feats <- cbind(feats,tot=apply(feats,1,sum))

outN <- rpart.forest.predict(t(datN@signal),hatN,augmented=T)

boxplot(outN[,"pred"]~surv[,2],xlab="Status",ylab="Expected hazard",col=c("dark green","dark red"))
x11(); boxplot(predict(hatN$forest[[1]])~surv[,2],xlab="Status",ylab="Expected hazard",col=c("dark green","dark red"))

hatN.2 <- rpart.forest.estimate( t(datN@signal), cls=surv, ntree=100 )
feats.2 <- hatN.2$features[apply(hatN.2$features,1,sum)>0,]
feats.2 <- cbind(feats.2,tot=apply(feats.2,1,sum))
outN.2 <- rpart.forest.predict(t(datN@signal),hatN.2,augmented=T)

boxplot(outN.2[,"pred"]~surv[,2],xlab="Status",ylab="Expected hazard",col=c("dark green","dark red"))

set.seed(123)
folds <- xval.select(ncol(datN),cls=surv[,2],nfolds=3)

DATtrn <- subset.resdata(datN,exptnames=exptnames(datN)[folds!=3])
CLStrn <- surv[folds!=3,]
DATtst <- subset.resdata(datN,exptnames=exptnames(datN)[folds==3])
CLStst <- surv[folds==3,]
HATtrn <- rpart.forest.estimate( t(DATtrn@signal), cls=CLStrn, ntree=100 )
OUTtrn <- rpart.forest.predict(t(DATtrn@signal), model=HATtrn, augmented=T)
boxplot(OUTtrn[,"pred"]~CLStrn[,2],xlab="Status",ylab="Expected hazard",col=c("dark green","dark red"))

OUTtst <- rpart.forest.predict(t(DATtst@signal), model=HATtrn, augmented=T)
boxplot(OUTtst[,"pred"]~CLStst[,2],xlab="Status",ylab="Expected hazard",col=c("dark green","dark red"))
