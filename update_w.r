library("Matrix")
#library("matrixcalc")
source("nrm.r")

QUIETW <- T
MAX_ITER <- 100
ABSTOL <- 1e-4
RELTOL <- 1e-2
RELAXPAR <- 1.2

update_w <- function(nsamples,A,b,Atb,AtA,Lapl,S,iS,nS,iM,iMS,L,U,z_init,u_init,rho,lam1,lam2) {
	
	n <- nrow(A)
	p <- ncol(A)
	
	x <- numeric(p)
	z <- numeric(p)
	u <- numeric(p)
	# Warm start z and u
	if(!is.null(z_init)) 
		z <- z_init
	if(!is.null(u_init)) 
		u <- u_init

	###ADMM solver
	
	# Cache the factorization
	#Z = t(A); Z = Z / sqrt(nsamples)
	#iMZ = backsolve(U, forwardsolve(L, Z))
	#D = diag(n) +  t(iMZ)%*%Z
	S = sqrt(lam2)*S

	history2 <- list()
	history2$objval <- list()
	history2$r_norm <- list()
	history2$s_norm <- list()
	history2$eps_pri <- list()
	history2$eps_dual <- list()
	
	#cat("\n")
	#print(paste("[update_w]norm(z)=",norm(z),", norm(u)=",norm(u),", norm(Atb)=",norm(Atb),sep=""))
	
	for (k2 in 1:MAX_ITER) {
		
		## x update
    		q = Atb + rho*(z - u) 
		##q = Atb + rho * S %*% ( t(S) %*% z - u)
		#nu = solve(D, t(iMZ)%*%q)
		#x = as.vector(backsolve(U, forwardsolve(L, q - Z%*%nu)))	
		nu = backsolve(U, forwardsolve(L, t(iMS)%*%q))
		x = as.vector(iM%*%(q - S%*%nu))
				
		## z update
    		zold = z
    		x_hat = RELAXPAR*x + (1 - RELAXPAR)*zold
		z = shrinkage(x_hat + u, lam1 / rho)
		#z <- shrinkage(x_hat + iS %*% u, lam1 / (rho * nS) )
		#FIXME
		#print(paste("[update_w]min(z)=",min(x_hat + u),", max(z)=",max(x_hat + u),", lam1=",lam1,sep=""))

		## u update
    		u = u + (x_hat - z)
    		#u <- u + t(S) %*% (x_hat - z)
		
		# Diagnostics, reporting, termination checks
		#history2$objval[[k2]] <- obj_func(A,b,Lapl,lam1,lam2,x,z)
		history2$r_norm[[k2]] <- norm(x-z)
		history2$s_norm[[k2]] <- norm(-rho*(z-zold))
		history2$eps_pri[[k2]] <- sqrt(p) * ABSTOL + RELTOL * max(norm(x), norm(-z)) + 0.0001
		history2$eps_dual[[k2]] <- sqrt(p) * ABSTOL + RELTOL * norm(rho*u) + 0.0001
		
		# rho update 
#		if (history2$r_norm[[k2]]>10L*history2$s_norm[[k2]]) {
#			rho_new = 2L*rho
#		} else if (history2$s_norm[[k2]]>10L*history2$r_norm[[k2]]) {
#			rho_new = rho / 2L
#		} else { 
#			rho_new = rho
#		}
#		scale = rho / rho_new
#		u = scale * u
#		rho = rho_new
#
		if (QUIETW==FALSE) 
			#print(paste("[update_w]k2=", ", objval=",history2$objval[[k2]], k2,", r_norm=", history2$r_norm[[k2]], ", eps_pri=", history2$eps_pri[[k2]], ", s_norm=", history2$s_norm[[k2]], ", eps_dual=", history2$eps_dual[[k2]], sep=""))
			 print(paste("[update_w]k2=",k2,", lam1=",lam1,", max(x)=",max(x),", norm(x)=",norm(x),", norm(z)=",norm(z),", norm(u)=",norm(u),sep=""))

		if (history2$r_norm[[k2]] <= history2$eps_pri[[k2]] && history2$s_norm[[k2]] <= history2$eps_dual[[k2]]) {
			#print(paste("[update_w]k2=", k2,", objval=", history2$objval[[k2]],", r_norm=", history2$r_norm[[k2]], ", eps_pri=", history2$eps_pri[[k2]], ", s_norm=", history2$s_norm[[k2]], ", eps_dual=", history2$eps_dual[[k2]], sep=""))
    			break
		}	
	}	
	#DEBUG
	#print(paste("[update_w]k2=",k2,", norm(z)=", norm(z), sep=""))
	
	obj = obj_func(A,b,Lapl,lam1,lam2,z,z)
	return(list(x=x,z=z,u=u,obj=obj))

}

obj_func = function(A, b, Lapl, lam1, lam2, x, z) {
   obj =  (1/2)*norm(b - A%*%x)^2 + (lam2/2)*t(x)%*%Lapl%*%x  + lam1*sum(abs(z))
   return(obj)		
}

shrinkage <- function(x, kappa) {
    #res <- pmax( 0, as.vector(x) - kappa ) - pmax( 0, -as.vector(x) - kappa )
    res <- pmax( 0, x - kappa ) - pmax( 0, -x - kappa )
    return(res)
}

factor_w <- function(M) {
    U <- chol( M )
    return(list(L=t(U),U=U))
}


