library("multicore")
source("admm_lasso.r")

QUIET <- T

admm_lasso_multi <- function(A,b,init_x,decomps,lam1,rho,relaxpar,sselection,bsize,ncores) {
	
	n <- dim(A)[1]
	p <- dim(A)[2]
	nblocks <- p/bsize
	print(paste("nblocks=", nblocks,sep=""))

	###ADMM solver
	
	history <- list()
	history$objval <- list()
	history$r_norm <- list()
	history$s_norm <- list()
	history$eps_pri <- list()
	history$eps_dual <- list()
	
	res_x <- vector("list", nblocks)
	
	#initializations
	ax_bar <- Matrix(numeric(n))
	z <- Matrix(0, nrow=n, ncol=nblocks)
	u <- Matrix(numeric(n))

	#warm start x
	for (j in 1:nblocks) {
		res_x[[j]]$x <- init_x[((j-1)*bsize+1):(j*bsize)]
		res_x[[j]]$ax <- A[, ((j-1)*bsize+1):(j*bsize)] %*% res_x[[j]]$x
		ax_bar <- ax_bar + res_x[[j]]$ax
    	}
	ax_bar <- ax_bar / nblocks
			
	z_bar <- (b + rho*(ax_bar + u)) / (nblocks+rho)
	
	for (j in 1:nblocks) 
		z[ , j] <- z_bar + res_x[[j]]$ax - ax_bar	
		
    	u <- u + ax_bar - z_bar
	
	for (k in 1:MAX_ITER) {

		#DEBUG
		#print(paste("[lasso]ADMM iteration ",k,sep=""))
		
		#x update in parallel
		cat("\n")
		res_x <- mclapply(1:nblocks, function(j) admm_lasso(A[ ,((j-1)*bsize+1):(j*bsize)], res_x[[j]]$ax + z_bar - ax_bar - u, res_x[[j]]$x, decomps[[j]]$L, decomps[[j]]$U, lam1, rho, relaxpar, sselection), mc.preschedule = FALSE, mc.cores = ncores)
		#print("length(res_x)")
		#print(length(res_x))
		#for (j in 1:nblocks) {
		#	print(paste("j=",j,sep=""))
		#	print(str(res_x[[j]]))
		#}
		aux <- res_x[[1]]$x
		if (nblocks>1) 
			for (j in 2:nblocks)
				aux <- rbind2(aux, res_x[[j]]$x)
		print(paste("[admm_lasso_multi]length(x)=", length(aux), ", norm(x)=", sqrt(sum(aux^2)),sep=""))
						
		z_bar_old <- z_bar
		ax_bar <- Matrix(numeric(n))
		for (j in 1:nblocks) 
			ax_bar <- ax_bar + res_x[[j]]$ax
		ax_bar <- ax_bar / nblocks
		ax_bar_hat <- relaxpar*ax_bar + (1-relaxpar)*z_bar_old
		
		#z update 
		z_bar <- (b + rho*(ax_bar_hat + u)) / (nblocks+rho)
		
		z_old <- z
		for (j in 1:nblocks) 
			z[ , j] <- z_bar + res_x[[j]]$ax - ax_bar #ax_bar_hat	
		
		#u update
    		u <- u + ax_bar_hat - z_bar
    		
		#DEBUG
		#print(paste("[lasso]norm(z_bar)=", sqrt(sum(z_bar^2)), ", norm(u)=", sqrt(sum(u^2)), sep=""))	
		
		#check stopping conditions
		p <- 0; s <- 0; ep1 <- 0; ep2 <- 0; es <- 0
		for (j in 1:nblocks) {
			ax_hat <- res_x[[j]]$ax
			p <- p + sum((ax_hat - z[ ,j])^2)

			Aj <- A[, ((j-1)*bsize+1):(j*bsize)]
			s <- s + sum((-rho*crossprod(Aj, (z[ ,j]-z_old[ ,j])))^2)

			ep1 <- ep1 + sum(ax_hat^2)
			ep2 <- ep2 + sum(z[ ,j]^2)

			es <- es + sum(crossprod(Aj, rho*u)^2)
		}

		#diagnostics, reporting, termination checks
		history$objval[[k]] <- opt_obj_multi(A,b,lam1,res_x,z_bar,nblocks)
		#history$r_norm[[k]] <- sqrt(p)
		history$r_norm[[k]] <- sqrt(nblocks) * sqrt(sum((z_bar - ax_bar)^2))
		history$s_norm[[k]] <- sqrt(s)
		history$eps_pri[[k]] <- sqrt(n) * ABSTOL + RELTOL * max(sqrt(ep1), sqrt(ep2))
		history$eps_dual[[k]] <- sqrt(bsize) * ABSTOL + RELTOL * sqrt(es)
		
		if (QUIET==FALSE)  
			print(paste("[admm_lasso_multi]k=", k,", r_norm=", history$r_norm[[k]], ", eps_pri=", history$eps_pri[[k]], ", s_norm=", history$s_norm[[k]], ", eps_dual=", history$eps_dual[[k]], ", objval=", history$objval[[k]], sep=""))

		if (history$r_norm[[k]]<history$eps_pri[[k]]&&history$s_norm[[k]]<history$eps_dual[[k]])
        		break
			 
	}#end for
	print(paste("[admm_lasso_multi]Num of ADMM iters ", k, sep=""))
		
	x <- res_x[[1]]$x
	if (nblocks>1) 
		for (j in 2:nblocks)
			x <- rbind2(x, res_x[[j]]$x)
		
	return(x)

}#end lasso

opt_obj_multi <- function(A, b, lambda, x, z, nblocks) {
   
   s_x <- 0
   #s_z <- Matrix(numeric(n))
   for (j in 1:nblocks) {
   	s_x <- s_x + sum(abs(x[[j]]$x))
	#s_z <- s_z + z[, j]
   }	
   #obj <-  (1/2)*sum((s_z-b)^2) + lambda*s_x
   obj <-  (1/2)*sum((nblocks*z-b)^2) + lambda*s_x
   
   return(obj)		
}


