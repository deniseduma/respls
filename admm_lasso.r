library("Matrix")

QUIET <- T
MAX_ITER <- 1000
ABSTOL <- 1e-3
RELTOL <- 1e-3

mu <- 10
pi_up_down <- 2

admm_lasso <- function(A,b,init_x,L,U,lam1,rho,relaxpar,sselection) {
	
	n <- dim(A)[1]
	p <- dim(A)[2]
	#cat("\t")
	#print("[admm_lasso]dim(A)")	
	#cat("\t")
	#print(dim(A))	
	
	###ADMM solver
	#save a matrix-vector multiply
	Atb <- crossprod(A, b)

	x <- Matrix(numeric(p))
	#z <- Matrix(numeric(p))
	z <- init_x
	u <- Matrix(numeric(p))

	history2 <- list()
	history2$objval <- list()
	history2$r_norm <- list()
	history2$s_norm <- list()
	history2$eps_pri <- list()
	history2$eps_dual <- list()
	
	for (k2 in 1:MAX_ITER) {

		#x update
    		q <- Atb + rho * (z - u)    # temporary value
    		if( n >= p ) {   # if skinny
			x <- backsolve(U,  forwardsolve(L, q))
    		} else {	# if fat
			x <- q/rho - crossprod(A, backsolve(U, forwardsolve(L, (A%*%q))))/(rho^2)
		}
		
		#DEBUG
		#print(paste("[",k,"]", "max(x)=", max(x), ", min(x)=", min(x), sep=""))
		#print("shrinkage")
		#print(lambda/rho)

		#z update
    		zold <- z
    		x_hat <- relaxpar*x + (1 - relaxpar)*zold
    		z <- shrinkage(x_hat + u, lam1/rho)

		#u update
    		u <- u + (x_hat - z)
    		
		#diagnostics, reporting, termination checks
		history2$objval[[k2]] <- opt_obj(A,b,lam1,x,z)
		history2$r_norm[[k2]] <- sqrt(sum((x-z)^2))
		history2$s_norm[[k2]] <- sqrt(sum((-rho*(z-zold))^2))
		history2$eps_pri[[k2]] <- sqrt(p) * ABSTOL + RELTOL * max(sqrt(sum(x^2)), sqrt(sum((-z)^2)))
		history2$eps_dual[[k2]] <- sqrt(p) * ABSTOL + RELTOL * sqrt(sum((rho*u)^2))
		
		if (QUIET==FALSE) {
			cat("\t")
			print(paste("[admm_lasso]k2=", k2,", r_norm=", history2$r_norm[[k2]], ", eps_pri=", history2$eps_pri[[k2]], ", s_norm=", history2$s_norm[[k2]], ", eps_dual=", history2$eps_dual[[k2]], sep=""))
		}
    		
		if (history2$r_norm[[k2]] < history2$eps_pri[[k2]] && history2$s_norm[[k2]] < history2$eps_dual[[k2]])
    			break
	
		#update penalty parameter
		if (history2$r_norm[[k2]] > mu * history2$s_norm[[k2]]) {
			rho <- pi_up_down * rho
			u <- u / pi_up_down
			res <- factor(A, rho)
			L <- res$L 
			U <- res$U
		} else 	if (history2$s_norm[[k2]] > mu * history2$r_norm[[k2]]) {
			rho <- rho / pi_up_down
			u <- u * pi_up_down
			res <- factor(A, rho)
			L <- res$L 
			U <- res$U
		}	

	}#end for
	#cat("\t")
	#print(paste("[admm_lasso]Num of ADMM iters ", k2, sep=""))
	
	ax <- A %*% z
	
	return(list(x=z, ax=ax))

}#end lasso

opt_obj <- function(A, b, lambda, x, z) {
   obj <-  (1/2)*sum((b - A%*%x)^2) + lambda*sum(abs(z))
   return(obj)		
}

shrinkage <- function(x, kappa) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    #cat("\n")
    #print(paste("kappa=",kappa,sep=""))
    #print(paste("[Before]min(x)=",min(x),", max(x)=",max(x),", nnzero(x)=",nnzero(x),sep=""))
    res <- Matrix((pmax( 0, as.vector(x) - kappa ) - pmax( 0, -as.vector(x) - kappa )),nrow=n,ncol=p)
    #print(paste("[After]min(x)=",min(x),", max(x)=",max(x),", nnzero(x)=",nnzero(x),sep=""))
    return(res)
}

factor <- function(X, rho) {
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    #print("str(U)")
    if ( n >= p ) {#if skinny
	U <- chol( crossprod(X) + rho*Diagonal(p) )
	#res <- Cholesky( crossprod(X) + rho*Diagonal(p) )
	#L <- as(res, "Matrix")	
    	#P <- as(res, "pMatrix")
    	#L <- t(P) %*% L
    } else {#if fat
	U <- chol( Diagonal(n) + (1/rho)*tcrossprod(X) )
    }
    #print(str(t(U)))

    return(list(L=t(U),U=U))
}
