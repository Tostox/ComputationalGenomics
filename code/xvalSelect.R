many2one.cls <- function( cls )
{
  # takes a matrix with one variable/vector per column and turns it
  # into a single variable/vector with as many values as the Cartesian
  # product of the column variables
  #
  if ( is.null(dim(cls)) || ncol(cls)==1 ) {
    cls2 <- as.numeric(match(drop(cls),sort(unique(drop(cls)))))-1
    levels(cls2) <- levels(cls)
    return(cls2)
  }
  # ELSE ..
  #
  levs <- sort(apply(unique(cls),1,paste,collapse="."))
  cls.new <- as.numeric( match( apply(cls,1,paste,collapse="."), levs) )-1
  levels(cls.new) <- levs
  cls.new
}
xvalSelect <- function( sample.size, nfolds, cls=NULL )
{
  # returns a sample.size-length vector of fold assignments (labels
  # btw 1 and nfolds). If cls is non-NULL, it carries out stratified
  # cross-validation. If you want to use more than one variable for
  # stratification, use many2one to convert several (discrete) variables
  # into one, which will be passed to cls
  #
  if ( nfolds<2 )
    stop( "nfolds should be at least 2" )
  if ( nfolds>sample.size )
    stop( "nfolds cannot be greater than sample size" )
  if ( nfolds<sample.size & sample.size%/%nfolds<2 )
    warning( nfolds-sample.size%%nfolds, " fold(s) w/ a single sample" )

  if ( nfolds==sample.size )
    return( 1:sample.size )
  
  if ( is.null(cls) )
  {
    fsize <- sample.size%/%nfolds
    rsize <- sample.size%%nfolds
    folds <- sample( c(rep(1:nfolds,fsize), if(rsize>0) sample(nfolds,rsize)) )
    return(folds)
  }
  else
  {
    if ( sample.size!=length(cls) )
      stop( "cls must be same length as sample size" )
    
    nstrata <- length( strata <- as.vector(sort(unique(cls))) )

    folds <- rep( NA, sample.size )
    for ( s in 1:nstrata )
    {
      idx <- cls==strata[s]
      if (sum(idx)==1) {
        folds[idx] <- sample(nfolds,1)
        next
      }
      folds[(1:sample.size)[idx]] <- xvalSelect( sum(idx), nfolds )
    }
    return(folds)
  }
}
