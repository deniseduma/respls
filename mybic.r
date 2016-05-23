
select_best_lambda <- function(ws, W, X0, X, Y, tvar) {
	
	n <- nrow(X)
	p <- ncol(X)
	if (n==1) {
		n <- n + 0.001
	}

	idx_lam <- 0
	best_w <- NULL
	min_bic <- 1e6
	for (j in 1:ncol(ws)) {
		
		w <- as(ws[, j], "sparseMatrix")
		nnz <- nnzero(w)
		if (nnz==0)
			next
								
		# Degree of freedom for BIC
		df <- nnz
		# -2log_likelihood
		err <- log(sum((Y - X%*%w)^2)/n) 
		# BIC
		crt_bic <- err + (log(n)/n) * df 
				
		if (crt_bic<min_bic) {
			idx_lam <- j
			best_w <- w 
			min_bic <- crt_bic
		}
		
		#DEBUG
		#if (j<=20) {
			print(paste("[bic]j=",j,", df=",df,", bic=",crt_bic,sep=""))
			#print(which(ws[,j]!=0))
			#print(ws[which(ws[,j]!=0)])
		#}		
	
	}#end for j
	
	BB <- W
	if (is.null(W)) {
		BB <- best_w
	} else {	
		BB <- cbind2(BB, best_w)
	}
	pW <- BB %*% solve( crossprod(BB), t(BB) )
	cpev <- sum( diag( crossprod(X0%*%pW) ) ) / tvar

	#DEBUG
	print(paste("[bic]tvar=", tvar,", cpev=", cpev,sep=""))
	print(paste("[bic]idx_lam=",idx_lam,", min_nnz=",nnzero(best_w),", min_bic=",min_bic,sep=""))
	print("[bic]nnzero indices")
	print(which(best_w!=0))

	return(list(idx_lam=idx_lam, bic=min_bic, cpev=cpev))
}
