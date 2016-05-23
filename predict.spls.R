
# return fit or coefficients

"predict.spls" <-
function( object, newx, type = c("fit","coefficient"), ... )
{
    # newx: matrix of predictors
    # type = "fit" or "coefficient"
    
    type <- match.arg(type)    
    betahat <- object$betahat
    #print("betahat")
    #print(betahat)
    #print("norm(beta1)")
    #print(norm(betahat[ ,1]))
    #print("norm(beta1)")
    #print(norm(betahat[ ,2]))
    x <- object$x
    A <- object$A
    p <- ncol(x)
    
    if ( type=="fit" )
    {        
        if ( missing(newx) )
        {
            pred <- x %*% betahat + matrix(1,nrow(x),1) %*% object$mu
        } else
        {
            if ( ncol(newx)!=p & ncol(newx)!=length(A) )
            { stop("The dimension of test dataset is inapproprite!") }
            if ( ncol(newx)==p ) { newx <- newx[,A,drop=FALSE] }
            newx <- scale( newx, object$meanx[A], object$normx[A] )
            pred <- newx %*% betahat[A,,drop=FALSE] + matrix(1,nrow(newx),1) %*% object$mu
        }
    }
    if ( type=="coefficient" ) { pred <- betahat }
    
    invisible(pred)
}

"coef.spls" <-
function( object, ... )
{
    predict.spls( object, type="coefficient", ... )
}
