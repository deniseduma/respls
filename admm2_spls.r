library("rdetools")

source("admm_lasso2_spls.r")
source("sampleTPoints.r")

admm2_spls <- function(X, Y, XtX, Lapl, x_init=NULL, nlambdas) {
		
	n <- nrow(X)
	p <- ncol(X)
	
	rho <- 1
	relaxpar <- 1.8
	
	x <- Matrix(numeric(p))
	u <- Matrix(numeric(p))
	if(!is.null(x_init)) 
		x <- x_init
	
	max_lam <- max(abs(crossprod(X, Y)))
	#if (n>p){
	#	min_lam <- 0.0001*max_lam
	#} else {
	#	min_lam <- 0.01*max_lam
	#}
	lambdas <- logspace(log10(0.01*max_lam),log10(0.95*max_lam),nlambdas)
	lambdas <- sort(abs(lambdas),decreasing=TRUE)
	
	print("[admm]lambdas")
	print(lambdas[c(1, length(lambdas))])
	
	xxs <- Matrix(0, nrow=p, ncol=length(lambdas))
	
	for(l in 1:length(lambdas)){
		#cat("\n")
		#print(paste('[', l, ']lambda=', lambdas[l],sep=""))
		
		lam1 <- lambdas[l] #0.8*lambdas[l]
		lam2 <- 0.4*lambdas[l]
		
		#cache the factorization
		res <- factor2(XtX,Lapl,lam2,rho)
		L <- res$L 
		U <- res$U
		
		res <- admm_lasso2_spls(X,Y,XtX,Lapl,L,U,x,u,lam1,lam2,rho,relaxpar)
		#temp <- (1 + lam2) * temp

		xxs[ ,l] <- res$x
				
		#warm start
		x <- res$x
		u <- res$u
	}#end for
		
	#plot regularization path
	#dev.new()
	#matplot(log10(lambdas), t(xxs), cex=1, pch=16, type="l", xaxt="n", xlim=c(log10(min_lam1), log10(max_lam1)), xlab ="lambda", ylab="Coeffs/Selection prob", main="Regularization path", cex.lab=1.25, cex.axis=1.25)
		
	return(xxs)

}	

