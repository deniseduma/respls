## Select best w using BIC
	cat("\n"); print("BIC")
	betas = list(); prs = list(); 
	idx_bic = 1; min_bic = 1e6; all_zero = 1
	for (j in 1:nlam) { 
		w = ws[ ,j]; idx_w = which(w!=0); nnz_w = length(idx_w)
		if (nnz_w==0) 
			next
		
		all_zero = 0 
			
		## OLS re-estimation of w
		#y = ys[ ,j]; sX = X[ ,idx_w,drop=FALSE]
		#sw = solve ( crossprod(sX) + 0.001*diag(nnz_w), t(sX)%*%y); w2[idx_w] = sw
		#err = (norm(y - X%*%w2)^2) / (q)
		
		## PLS re-estimation of w 
        	idx = union( idx_W, idx_w )
        	sX0 = X0[,idx,drop=FALSE]
        	plsfit = pls::plsr( Y0~sX0, ncomp=min(k,length(idx)), method="simpls", scale=F)
        	beta = matrix( 0, p, q ); beta[idx,] = matrix( coef(plsfit), length(idx), q )
        	betas[[j]] = beta; prs[[j]] = plsfit$projection
		err = (norm(Y0-X0%*%beta)^2) / (n*q)
		
		## df estimation
		## Method 1
		dfterm = (log(n*q) / (n*q)) * nnz_w
		## Method 2
        	#sX = X[,idx_w,drop=FALSE]
		#lam = lambdas[j]; lam2 = alpha2*lam
		#DF = sX %*% solve( t(sX)%*%sX + lam2*L[idx_w, idx_w] + lam2*diag(nnz_w), t(sX))
		#df =  sum(diag(DF)); dfterm = (log(n*q) / (n*q))*df
			
		#plty = 2*0.5*log(choose(p, nnz_w))
		bic = log(err) + dfterm # + plty
		print(paste("j=",j,", nnz=",nnz_w,", err=",err,", log(err)=",log(err),", dfterm=",dfterm, ", bic=",bic,sep=""))

		if (bic<min_bic) {
			idx_bic = j; min_bic = bic
		}
	}
		
	if (all_zero) {
		print(paste("[opt]No non-zero found at PC ",k,sep=""))
		break
	}	

	w = ws[ ,idx_bic]
	v = vs[ ,idx_bic]

