## function: PLOT DISTN
##
plot.distn <- function(qfun, dfun, n, type="l", pch=20, lty=1, add=FALSE, xlab=NULL, ylab=NULL,
                       xlim=NULL,ylim=NULL,main=NULL, col="black", silent=TRUE, ...)
{
  p <- sort(c(0.000001,(1:(n-1))/n,.999999))
  q <- qfun(p,...)
  
  if ( add ) {
    if (type=="p")
      points( q, dfun(q,...), pch=pch, col=col )
    else
      lines( q, dfun(q,...), lty=lty, col=col )
  }
  else
      plot(q, dfun(q,...), type=type, lty=lty, pch=pch, col=col,
           xlab=xlab, ylab=ylab, main=main, xlim=xlim, ylim=ylim )

  ## this is useful to determine plot boundaries
  if (!silent) return ( list(xlim=c(q[1],q[length(q)]),ylim=c(0,max(dfun(q,...)))) )
}

## function: PLOT GAMMA
##
plot.gamma <- function( n, shape, scale=1, add=F, lty=1, ... )
{
  p <- (1:(n-1))/n
  q <- qgamma(p,shape=shape,scale=scale)
  if ( add )
    lines( q, dgamma(q,shape=shape,scale=scale), lty=lty, ... )
  else
    plot( q, dgamma(q,shape=shape,scale=scale), lty=lty, ... )
}
## function: PLOT LNORM
##
plot.lnorm <- function( n, meanlog, sdlog, add=F, lty=1, ... )
{
  p <- (1:(n-1))/n
  q <- qlnorm(p,meanlog=meanlog,sdlog=sdlog)
  if ( add )
    lines( q, dlnorm(q,meanlog=meanlog,sdlog=sdlog), lty=lty, ... )
  else
    plot( q, dlnorm(q,meanlog=meanlog,sdlog=sdlog), type="l", lty=lty, ... )
}
# TO PLOT A DISTRIBUTION
#
# p <- (1:100)/100    # in general, (1:n-1)/n
# lines( p, d<probfun>(q<probfun>)(p,...))
#
# for example:

if ( FALSE )
{
    p <- (1:999)/1000 
    shape <- 20.85198;scale <- 1/10.54734;plot( qgamma(p,shape=shape,scale=scale), dgamma( qgamma(p,shape=shape,scale=scale),shape=shape,scale=scale ),type="l", lty="solid" )
    
    plot( qnorm(p,mean=0,sd=1), dnorm( qnorm(p,mean=0,sd=1),mean=0,sd=1 ),type="l", lty="solid" )
    
    
    ## to multiply a vector my a matrix's rows
    ##
    t( vector * t(matrix) )

    ## multiple plots in a page
    ##
    par(mfrow=c(rows,cols))
    par(mfcol=c(cols,rows))
    
    ## histogram plus empirical density
    ##
    hist(x,probability=T,ylim=c(0,max(density(x)$y))); lines(density(x))
}
