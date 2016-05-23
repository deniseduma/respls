library("Matrix")
#library("matrixcalc")

QUIET <- T
MAX_ITER <- 1000
ABSTOL <- 1e-4
RELTOL <- 1e-2

mu <- 10

respls_admm_lasso <- function(A,b,Atb,Lapl,Le,M,z_init=NULL,u_init=NULL,lam,alpha,relaxpar) {
	
	q <- nrow(A)
	p <- ncol(A)
	
	rho <- 1 
	lam1 <- alpha*lam
	lam2 <- (1-alpha)*lam
	
	# Pre-computations 
	#DZM <- init_matrices(M,t(A),Diagonal(q))
	Z <- t(A)
	# V2
	#M <- Le$vectors%*%Diagonal(p, 1/(2*lam2*Le$values+rho*1))%*%t(Le$vectors)
	#Mz <- M%*%Z
	#D <- Diagonal(q) + t(Mz)%*%Z
	
	# V1
	M <- 2*lam2*Lapl + rho*Diagonal(p)
	Mz <- solve(M, Z)
	D <- Diagonal(q) + t(Mz)%*%Z
	
	###ADMM solver
	
	# Save a matrix-vector multiply
	#Atb <- crossprod(A, b)
	
	x <- Matrix(numeric(p))
	# warm start
	z <- z_init
	u <- u_init
	
	history2 <- list()
	history2$objval <- list()
	history2$r_norm <- list()
	history2$s_norm <- list()
	history2$eps_pri <- list()
	history2$eps_dual <- list()
	
	for (k2 in 1:MAX_ITER) {
		
		# x update
    		q <- Atb + rho * (z - u) 
		#nu <- DZM%*%q
		nu <- solve(D, t(Mz)%*%q)
		x <- solve(M, q - Z%*%nu)	
		#x <- M%*%(q - Z%*%nu)	

		# z update
    		zold <- z
    		x_hat <- relaxpar*x + (1 - relaxpar)*zold
    		z <- shrinkage(x_hat + u, lam1/rho)
		
		# u update
    		u <- u + (x_hat - z)
				    		
		#diagnostics, reporting, termination checks
		history2$objval[[k2]] <- opt_obj(A,b,Lapl,lam1,lam2,x,z)
		history2$r_norm[[k2]] <- sqrt(sum((x-z)^2))
		history2$s_norm[[k2]] <- sqrt(sum((-rho*(z-zold))^2))
		history2$eps_pri[[k2]] <- sqrt(p) * ABSTOL + RELTOL * max(sqrt(sum(x^2)), sqrt(sum((-z)^2)))
		history2$eps_dual[[k2]] <- sqrt(p) * ABSTOL + RELTOL * sqrt(sum((rho*u)^2))
		
		if (QUIET==FALSE) 
			print(paste("[respls_admm_lasso]k2=", k2,", r_norm=", history2$r_norm[[k2]], ", eps_pri=", history2$eps_pri[[k2]], ", s_norm=", history2$s_norm[[k2]], ", eps_dual=", history2$eps_dual[[k2]], sep=""))
    		
		if (history2$r_norm[[k2]] < history2$eps_pri[[k2]] && history2$s_norm[[k2]] < history2$eps_dual[[k2]])
    			break
	
		#update penalty parameter
#		if (history2$r_norm[[k2]] >= mu * history2$s_norm[[k2]]) {
#			rho <- pi_up_down * rho
#			u <- u / pi_up_down
#			res <- factor2(AtA,Lapl,lam2,rho)
#			L <- res$L 
#			U <- res$U
			#DEBUG
#			print(paste("[respls_admm_lasso]r_norm is bigger ",history2$r_norm[[k2]], ", vs ", history2$s_norm[[k2]],sep="") )
#		} else 	if (history2$s_norm[[k2]] > mu * history2$r_norm[[k2]]) {
#			rho <- rho / pi_up_down
#			u <- u * pi_up_down
#			res <- factor2(AtA,Lapl,lam2,rho)
#			L <- res$L 
#			U <- res$U
			#DEBUG
#			print(paste("[respls_admm_lasso]s_norm is bigger ",history2$s_norm[[k2]], ", vs ", history2$r_norm[[k2]],sep="") )
#		}	

	}
	
	return(list(z=z, u=u))

}

opt_obj <- function(A, b, Lapl, lam1, lam2, x, z) {
   obj <-  (1/2)*sum((b - A%*%x)^2) + lam2*t(x)%*%Lapl%*%x  + lam1*sum(abs(z))
   return(obj)		
}

shrinkage <- function(x, kappa) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    res <- Matrix((pmax( 0, as.vector(x) - kappa ) - pmax( 0, -as.vector(x) - kappa )),nrow=n,ncol=p)
    return(res)
}

#factor_w <- function(K,Lapl,lam2,rho) {
#    U <- chol( K + 2*lam2*Lapl + rho*Diagonal(p) )
#    return(list(L=t(U),U=U))
#}


