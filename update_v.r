#library("rdetools")
library("matrixcalc")

source("nrm.r")

QUIETV <- T
MAX_ITER <- 100
ABSTOL <- 1e-4
RELTOL <- 1e-2
RELAXPAR <- 1.2

update_v <- function(A, b, Atb, AtA, V, iM, iMV, VtiMV, z_init,u_init,rho) {
	
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
	#A <- A / sqrt(n)

	history2 <- list()
	history2$objval <- list()
	history2$r_norm <- list()
	history2$s_norm <- list()
	history2$eps_pri <- list()
	history2$eps_dual <- list()
	
	#DEBUG
	#cat("\n")
	#print(paste("?[update_v]norm(z)=",norm(z),", norm(u)=",norm(u),", norm(Atb)=",norm(Atb),sep=""))
	
	for (k2 in 1:MAX_ITER) {
		
		# x update
    		q <- Atb + rho * (z - u) 
		if (is.null(V)) { #first PC
			## Method 1: Cholesky decomposition
			#x = as.vector(backsolve(U, forwardsolve(L, q)))
			## Method 2: Woodbury formula
			x = as.vector(iM%*%q)
			## Method 3: Woodbury formula better computation
			#myW  = solve(diag(n) + (1/rho)*A%*%t(A), A %*% q)
			#x = (1/rho)*q - (1/rho^2)*t(A)%*%myW
		} else { #remaining PCs
			## Method 1: Cholesky decomposition
			#nu = solve(D, t(Vhat) %*% q)
			#x = as.vector(backsolve(U, forwardsolve(L, q - V%*%nu)))	
			## Method 2: Woodbury formula
			nu = solve(VtiMV, t(iMV)%*%q)
			x = as.vector(iM%*%(q - V%*%nu))
			## Method 3: Woodbury formula better computation
			#nu = solve(D, t(Vhat)%*%q)
			#myq = q - V%*%nu 
			#myW  = solve(diag(n) + (1/rho)*A%*%t(A), A%*%myq)
			#x = (1/rho)*myq - (1/rho^2)*t(A)%*%myW
		}
		
		# z update
    		zold <- z
    		x_hat <- RELAXPAR*x + (1 - RELAXPAR)*zold
    		z <- x_hat + u
		#if (norm(z)>1)
			z <- z / norm(z)
				
		# u update
    		u <- u + (x_hat - z)
				
		# Diagnostics, reporting, termination checks
		#history2$objval[[k2]] <- obj_func_v(A,b,V,nu,x,z)
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
	
		if (QUIETV==FALSE) 
			print(paste("[update_v]k2=", k2,", r_norm=", history2$r_norm[[k2]], ", eps_pri=", history2$eps_pri[[k2]], ", s_norm=", history2$s_norm[[k2]], ", eps_dual=", history2$eps_dual[[k2]], sep=""))
    		
		if (history2$r_norm[[k2]] <= history2$eps_pri[[k2]] && history2$s_norm[[k2]] <= history2$eps_dual[[k2]]) {
			#print(paste("[update_v]k2=", k2,", r_norm=", history2$r_norm[[k2]], ", eps_pri=", history2$eps_pri[[k2]], ", s_norm=", history2$s_norm[[k2]], ", eps_dual=", history2$eps_dual[[k2]], sep=""))
    			break
		}	
	
	}
	#print(paste("[update_v]k2=",k2,", norm(z)=",norm(z),", norm(u)=",norm(u),", norm(x)=",norm(x),sep=""))
	#print(paste("[update_v]k2=",k2,sep=""))
	
	obj = obj_func_v(A,b,V,nu,z,z)
	return(list(x=x,z=z,u=u,obj=obj))
}	

obj_func_v = function(A, b, V, nu, x, z) {
   if (is.null(V))
   	obj =  (1/2)*norm(b - A%*%x)^2
   else
   	obj =  (1/2)*norm(b - A%*%x)^2 + t(nu)%*%(t(V)%*%x)
   return(obj)		
}


factor_v <- function(M) {
    U <- chol( M )
    return(list(L=t(U),U=U))
}


