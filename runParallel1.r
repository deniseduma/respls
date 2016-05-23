library(MASS)
library(gplots)
library(RColorBrewer)

source("ust.R")
source("spls.R")
source("correctp.R")
source("spls.dv.R")
source("predict.spls.R")

source("opt.r")

runParallel1 <- function(i,n,p,p1,p2,q,yprop,mu1,mu2,npcs,nlam,V0,L,outdir,method) {
	
	cat("\n")
	print(paste("Iteration ", i, sep=""))
		
	# Generate X
	X <- mvrnorm(n, rep(0, p), diag(p))
	newX <- mvrnorm(n, rep(0, p), diag(p)) 
	
	# Simulate first p1 genes
	X[, 1:p1] = X[, 1:p1] + mu1
	newX[, 1:p1] = newX[, 1:p1] + mu1
	# Simulate next p2 genes
	X[, (p1+1):(p1+p2)] = X[, (p1+1):(p1+p2)] + mu2
	newX[, (p1+1):(p1+p2)] = newX[, (p1+1):(p1+p2)] + mu2

	cnames = ("X1")
	for (j in 2:p)
		cnames = c(cnames, paste("X",j,sep=""))
	colnames(X) = cnames	

	# Generate y
	#if (q==1) { 
		Y <- yprop*mu1 + (1-yprop)*mu2 + rnorm(n,0,1)
		newY <- yprop*mu1 + (1-yprop)*mu2 + rnorm(n,0,1)
	#} else if (q==2) {
	#	y1 <- mu1 + rnorm(n,0,1)
	#	y2 <- mu2 + rnorm(n,0,1)
	#	newy1 <- mu1 + rnorm(n,0,1)
	#	newy2 <- mu2 + rnorm(n,0,1)
	#	Y <- cbind(y1, y2)
	#	newY <- cbind(newy1, newy2)
	#}

	##Measure execution time
	start <- proc.time()
	if (method == "spls") {
		res = spls(X, Y, K=q, eta=0.5, kappa=0.5, select="pls2", fit="simpls", scale.x=TRUE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=FALSE)
	} else if (method == "respls") {
		res = opt(X, Y, L, NULL, NULL, 0, q, nlam,outdir)
	}
	#Measure execution time
	elapsed <- proc.time() - start
	
	if (method == "spls") {
		#betapls <- res$betahat
		#W = sapply(res$ws, unlist)
		W = as.matrix(res$ws[[1]])
		if (length(res$ws)>1)
			for (j in 2:length(res$ws)) 
				W = cbind(W, res$ws[[j]])
		# Normalize direction vectors 
		W = W %*% diag( 1/sqrt(colSums(W^2)), nrow=ncol(W), ncol=ncol(W) ) 
	} else if (method == "respls") {
		W = res$W
	}

	# Prediction
	if (method == "spls") {
		predY = predict.spls( res, newX, type = "fit" )
	} else if (method == "respls") {
		predY = predict.respls( res, newX)
	}
	# Mean squared prediction error
    	one <- matrix(1,1,n)
	mspe <- drop (one %*% ((newY - predY)^2)) / n

	# TPR and FPR
	tp = matrix(0,npcs,nlam)
	fp = matrix(0,npcs,nlam)
	tn = matrix(0,npcs,nlam)
	fn = matrix(0,npcs,nlam)

	recall = matrix(0, npcs, nlam)
	prec = matrix(0, npcs, nlam)

	tp2 = numeric(npcs)
	fp2 = numeric(npcs)
	tn2 = numeric(npcs)
	fn2 = numeric(npcs)
	
	recall2 = numeric(npcs)
	prec2 = numeric(npcs)
	
	if (method == "respls") {
		var = res$var
		cpev = res$cpev
	}
		
	cat("\n")
	best_match = numeric(npcs)
	V_hat2 = matrix(0,p,npcs)
	V_hat = matrix(0,p,(npcs*nlam))
	for (j in 1:npcs) {
		min_j = 1
		min_ang = 1e6
		for (j2 in 1:ncol(W)) { 
			ang = acos(abs(crossprod(V0[, j, drop=FALSE], W[ ,j2, drop=FALSE])))
			if (ang<min_ang) { 
				min_ang = ang
				min_j = j2
			}
		}
			
		#DEBUG
		print(paste("[runParallel1]iter=",i,", j=",j,", min_j=",min_j,", min_ang=",min_ang,sep=""))
		
		best_match[j] = min_j
		V_hat2[, j] =  W[ ,min_j]
		if (method == "respls") {
			V_hat[, ((j-1)*nlam + 1):(j*nlam)] =  res$bW[ ,((min_j-1)*nlam + 1):(min_j*nlam)]
		}	
	}
	
	# Were components correctly split?
	csplit <- length(unique(best_match)) / length(best_match)

	if (method == "respls") {
		for (j in 1:npcs) {
			V_hat_crt <- V_hat[, ((j-1)*nlam+1):(j*nlam)]
			for (l in 1:nlam) {
				tp[j,l] = sum(V0[,j] & V_hat_crt[,l])
				fn[j,l] = sum(V0[,j]!=0) - tp[j,l]
				fp[j,l] = sum((!V0[,j]) & V_hat_crt[,l])
				tn[j,l] = sum(!V0[,j]) - fp[j,l]
				recall[j,l] = tp[j, l]/(tp[j, l]+fn[j, l])
				prec[j,l] = tp[j, l]/(tp[j, l]+fp[j, l])
			}
		}
	}

	for (j in 1:npcs) {
		tp2[j] = sum(V0[,j] & V_hat2[,j])
		fn2[j] = sum(V0[,j]!=0) - tp2[j]
		fp2[j] = sum((!V0[,j]) & V_hat2[,j])
		tn2[j] = sum(!V0[,j]) - fp2[j]
		recall2[j] = tp2[j]/(tp2[j]+fn2[j])
		prec2[j] = tp2[j]/(tp2[j]+fp2[j])
	}	
	
	# Total FPR
	v0 <- V0[,1]		
	if (npcs>1)
		for (j in 2:npcs) {
			v0 <- v0 | V0[,j]
		}
	#print("v0")
	#print(which(v0!=0))
	w <- W[,1]
	if (npcs>1)
		for (j in 2:ncol(W)) {
			w <- w | W[,j]
		}
	#print("w")
	#print(which(w!=0))
	t_tp = sum(v0 & w)
	t_fp = sum((!v0) & w)
	t_tn = sum(!v0) - t_fp
	t_prec = t_tp/(t_tp+t_fp)
	

	if (method == "spls") {
		res = list(X=res$X,Y=res$Y,betapls=res$betapls,elapsed=elapsed,mspe=mspe,recall2=recall2,prec2=prec2,t_prec=t_prec,csplit=csplit)
	} else if (method == "respls") {
		res <- list(X=res$X,Y=res$Y,P=res$P,Q=res$Q,betapls=res$betapls,elapsed=elapsed,var=var,cpev=cpev,mspe=mspe,recall=recall,prec=prec,recall2=recall2,prec2=prec2,t_prec=t_prec,csplit=csplit)
	}
	
	return(res)

}
