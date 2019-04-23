sample.size <- function(theta0, theta1, sigma, alpha=.05, beta=.10, two.sided=T )
{
  z.alpha <- qnorm(alpha/2,mean=0,sd=sigma,lower.tail=F)
  z.beta <- qnorm(beta/2,mean=0,sd=sigma,lower.tail=T)

  n <- (sigma * ( z.alpha - z.beta ) / ( theta1 - theta0 ) )^2
  return( ceiling(if (two.sided) 2*n else n ) )
}


SIGMA <- seq(0.5,1,.05)
DELTA <- seq(0.1,0.5,0.1)
SS <- matrix(NA,nrow=length(DELTA),ncol=length(SIGMA),
             dimnames=list(DELTA,SIGMA))
for ( i in 1:length(DELTA) ) {
  for ( j in 1:length(SIGMA) ) {
    SS[i,j] <- sample.size(theta0=0,theta1=DELTA[i],sigma=SIGMA[j])
  }
}
tmp <- matplot(SIGMA,t(SS),pch=".",type="l",lty=1,log="xy",ylab="sample size",ylim=c(6,2000))
legend("topleft",legend=paste("delta=",DELTA,sep=""),col=c("light blue","blue","green","red","black"),pch=20,lty=1)

SIGMA <- seq(0.25,1,.05)
DELTA <- seq(0.1,0.5,0.1)
SS <- matrix(NA,nrow=length(DELTA),ncol=length(SIGMA),dimnames=list(DELTA,SIGMA))
for ( i in 1:length(DELTA) ) {
  for ( j in 1:length(SIGMA) ) {
    #sigma <- DELTA[i]*SIGMA[j]
    sigma <- SIGMA[j]
    SS[i,j] <- ceiling(power.t.test(power=.9,delta=DELTA[i],sd=sigma,sig.level=0.05)$n)
  }
}
COL <- c("light blue","blue","green","red","black","orange","dark green")
matplot(SIGMA,t(SS),pch=".",col=COL,type="l",lty=1,log="xy",ylab="sample size",ylim=c(6,2000),xaxt=F)
legend("topleft",legend=paste("delta=",DELTA,sep=""),col=COL[1:length(DELTA)],pch=20,lty=1)
