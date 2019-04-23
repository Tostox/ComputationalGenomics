SUM <- read.table("/home/radon00/smonti/wi/rforests/Summary.transposed.short.txt",header=T,sep="\t")
SUM <- read.table("/home/radon00/smonti/wi/rforests/summary.20060801.xls",header=T,sep="\t")
DAT <- as.matrix(SUM[,-(1:2)])
rownames(DAT) <- SUM[,1]
COL <- c("red","green","blue","orange","purple","dark green","dark red","magenta")

pdf( "summary.pdf")
xi <- 1
yi <- seq(3,17,2)
yi <- seq(3,15,2)
for (i in 1:length(yi)) {
  lims <- c(0,max(DAT[xi,],DAT[yi[i],],na.rm=T))
  tscore <- t.test(DAT[xi,],DAT[yi[i],],paired=T,alternative="t")
  plot(DAT[xi,],DAT[yi[i],],pch=19,col=COL[i],xlab=rownames(DAT)[xi],ylab=rownames(DAT)[yi[i]],
       main=paste("ERR (t=",round(tscore$statistic,2),", p=",format.pval(tscore$p.value),")",sep=""),
       xlim=lims,ylim=lims)
  abline(0,1)
  #scan()
}
xi <- 2
yi <- seq(3,15,2)+1
for (i in 1:length(yi)) {
  cat(".")
  lims <- c(0,max(DAT[xi,],DAT[yi[i],],na.rm=T))
  tscore <- t.test(DAT[xi,],DAT[yi[i],],paired=T,alternative="t")
  plot(DAT[xi,],DAT[yi[i],],pch=19,col=COL[i],xlab=rownames(DAT)[xi-1],ylab=rownames(DAT)[yi[i]-1],
       main=paste("RERR (t=",round(tscore$statistic,2),", p=",format.pval(tscore$p.value),")",sep=""),
       xlim=lims,ylim=lims)
  abline(0,1)
  #scan()
}
dev.off()
