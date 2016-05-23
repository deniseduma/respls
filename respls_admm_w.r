library("rdetools")

source("respls_admm_lasso.r")
source("sampleTPoints.r")

alpha <- 0.8
RELAXPAR <- 1.2

respls_admm_w <- function(A, b, Lapl,Le,M,z_init=NULL, u_init=NULL, nlambdas) {
		
	q <- nrow(A)
	p <- ncol(A)
	
	z <- Matrix(numeric(p))
	u <- Matrix(numeric(p))
	if(!is.null(z_init)) 
		z <- z_init
	if(!is.null(u_init)) 
		u <- u_init
	
	# Generate lambda sequence
	Atb <- crossprod(A, b)
	max_lam <- (1/alpha)*max(abs(Atb))
	#if (n>p){
	#	min_lam <- 0.0001*max_lam
	#} else {
	#	min_lam <- 0.01*max_lam
	#}
	lambdas <- logspace(log10(0.001*max_lam),log10(0.9*max_lam),nlambdas)
	lambdas <- sort(abs(lambdas),decreasing=TRUE)
	
	print("[admm]lambdas")
	print(lambdas[c(1, length(lambdas))])
	
	ws <- Matrix(0, nrow=p, ncol=length(lambdas))
	us <- Matrix(0, nrow=p, ncol=length(lambdas))
	
	for(l in 1:length(lambdas)){
		
		lam <- lambdas[l]
		print(paste('[respls_admm_w]lam=',lam,sep=""))
		res <- respls_admm_lasso(A,b,Atb,Lapl,Le,M,z,u,lam,alpha,RELAXPAR)
		
		# Warm start
		z <- res$z
		u <- res$u

		ws[ ,l] <- z
		us[ ,l] <- u
	}
	return(list(ws=ws, us=us))
}	


