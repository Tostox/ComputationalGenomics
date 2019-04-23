set.seed(123)
N <- 500

X0 <- rnorm(N,-1,1)
Y0 <- rnorm(N,-1,1)
X1 <- rnorm(N,1,1)
Y1 <- rnorm(N,1,1)

xlim <- c(min(X0,X1)-1,max(X0,X1)+1)
ylim <- c(min(Y0,Y1)-1,max(Y0,Y1)+1)

pch <- 1
png("redundancy.red.png")
plot(X0, Y0, xlim=xlim, ylim=ylim, pch=pch, col="red",
     xlab="X",
     ylab="Y")
dev.off()
png("redundancy.green.png")
plot(X1, Y1, xlim=xlim, ylim=ylim, pch=pch, col="green",
     xlab="X",
     ylab="Y")
dev.off()

png("redundancy.bare.png")
plot(X0, Y0, xlim=xlim, ylim=ylim, pch=pch, col="red",
     xlab="X",
     ylab="Y")
points( X1, Y1, pch=pch, col="green")
dev.off()

png("redundancy1D.png")
pch <- 1
plot(X0, Y0, xlim=xlim, ylim=ylim, pch=pch, col="red",
     xlab=paste("mean(X|green)-mean(X|red)=",round(mean(X1)-mean(X0),2),sep=""),
     ylab=paste("mean(Y|green)-mean(Y|red)=",round(mean(Y1)-mean(Y0),2),sep=""))
points( X1, Y1, pch=pch, col="green")

points( X1, rep(ylim[1],N), pch="|", col="green")
points( X0, rep(ylim[1],N), pch="|", col="red")
points( c(mean(X0),mean(X1)), rep(ylim[1]+.3,2), pch=19, col="black" )
lines(  c(mean(X0),mean(X1)), rep(ylim[1]+.3,2), col="black", lwd=2 )
points( rep(xlim[1],N), Y1, pch="_", col="green" )
points( rep(xlim[1],N), Y0, pch="_", col="red" )
points( rep(xlim[1],2)+.2, c(mean(Y0),mean(Y1)), pch=19, col="black" )
lines(  rep(xlim[1],2)+.2, c(mean(Y0),mean(Y1)), col="black", lwd=2 )
dev.off()

png("redundancy.png")
pch <- 1
plot(X0, Y0, xlim=xlim, ylim=ylim, pch=pch, col="red",
     xlab=paste("mean(X|green)-mean(X|red)=",round(mean(X1)-mean(X0),2),sep=""),
     ylab=paste("mean(Y|green)-mean(Y|red)=",round(mean(Y1)-mean(Y0),2),sep=""))
points( X1, Y1, pch=pch, col="green")

points( X1, rep(ylim[1],N), pch="|", col="green")
points( X0, rep(ylim[1],N), pch="|", col="red")
points( c(mean(X0),mean(X1)), rep(ylim[1]+.3,2), pch=19, col="black" )
lines(  c(mean(X0),mean(X1)), rep(ylim[1]+.3,2), col="black", lwd=2 )

points( rep(xlim[1],N), Y1, pch="_", col="green" )
points( rep(xlim[1],N), Y0, pch="_", col="red" )
points( rep(xlim[1],2)+.2, c(mean(Y0),mean(Y1)), pch=19, col="black" )
lines(  rep(xlim[1],2)+.2, c(mean(Y0),mean(Y1)), col="black", lwd=2 )

points(mean(X0),mean(Y0),col="black")
points(mean(X0),mean(Y0),col="black",pch=19)
points(mean(X1),mean(Y1),col="black",pch=19)
lines(c(mean(X0),mean(X1)),c(mean(Y0),mean(Y1)),lwd=2)
dev.off()

rot <- function(theta) matrix( c(cos(theta),-sin(theta),sin(theta),cos(theta)), 2, 2)

dat <- rbind(cbind(X0,Y0),cbind(X1,Y1))

datR <- dat %*% rot(-pi/4)
dst <- abs(mean(datR[1:N,1])-mean(datR[(N+1):(2*N),1]))

png("redundancy45.png")
plot(datR[,1],datR[,2],pch=pch,xlab=paste("mean(X'|green)-mean(X'|red)=",round(dst,2),sep=""),ylab="Y'",xlim=xlim,ylim=ylim,main="45-degree rotation")
points(datR[1:N,1],datR[1:N,2],col="red",pch=pch)
points(datR[(N+1):(2*N),1],datR[(N+1):(2*N),2],col="green",pch=pch)
dev.off()

COR <- cbind(c("cor(X,Y|red)",
               "cor(X,Y|green)",
               "cor(X,Y)"),
             round(c(cor(X0,Y0),cor(X1,Y1),cor(c(X0,X1),c(Y0,Y1))),2))
my.write.matrix(COR,col.names=F,justify="r")
