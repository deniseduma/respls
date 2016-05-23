library(MASS)
library(gplots)
#library(RColorBrewer)

source("ust.R")
source("spls.R")
source("correctp.R")
source("spls.dv.R")
source("predict.spls.R")

source("opt.r")

runParallel3 <- function(i,n,p,q,npcs,nlam,C,V0,L,S,invS,nrmS,pcs,outdir,method) {
	
	cat("\n")
	print(paste("Iteration ", i, sep=""))
	
	#Generate train and test X
	X = mvrnorm(n, rep(0, p), C)
	newX = mvrnorm(n, rep(0, p), C)
	
	cnames = ("X1")
	for (j in 2:p)
		cnames = c(cnames, paste("X",j,sep=""))
	colnames(X) = cnames	

	# Generate Y
	
	## Method 1
	#U <- X %*% V0
	#newU <- newX %*% V0
	## Method 2
	r = svd(X); U = r$u[ ,1:npcs,drop=F]; V = r$v[ ,1:npcs,drop=F]; Sig = diag(r$d[1:npcs],npcs,npcs); Z = U %*% Sig 
	newr = svd(newX); newU = newr$u[ ,1:npcs,drop=F]; newV = newr$v[ ,1:npcs,drop=F]; newSig = diag(newr$d[1:npcs],npcs,npcs); newZ = newU %*% newSig
	print(paste("npcs ",npcs,sep="")); print("Sig"); print(Sig)
	
	cat("\n")
	best_match = numeric(npcs)
	for (j in 1:ncol(V0)) {
		min_j = 1; min_ang = 1e6
		for (j2 in 1:ncol(V)) { 
			ang = acos( abs ( crossprod( V0[ ,j], V[ ,j2] ) ) / ( norm(V0[ ,j])*norm(V[ ,j2]) ) ) 
			if (ang<min_ang) { 
				min_ang = ang
				min_j = j2				
			}
		}
		print(paste("[runParallel3][svd]j=", j, ", min_j=", min_j,", min_ang=",min_ang, sep=""))
		best_match[j] = min_j
	}
	#best_pcs = best_match[pcs]; print("best_pcs"); print(best_pcs)
	#U = as.matrix(U[ ,best_pcs]); Sig = as.matrix(Sig[best_pcs, best_pcs]); Z = U %*% Sig
	#newU = as.matrix(newU[ ,best_pcs]); newSig = as.matrix(newS[best_pcs, best_pcs]); newZ = newU %*% newSig
	cat("\n"); print("dim(U)"); print(dim(U)); print("dim(Sig)"); print(dim(Sig))
	
	#q1 = c(0.8, 0.04, 0.04, 0.04, 0.04) 
	#q2 = c(0.04, 0.8, 0.04, 0.04, 0.04) 
	#q3 = c(0.04, 0.pi 0.8, 0.04, 0.04) 
	#q4 = c(0.04, 0.04, 0.04, 0.8, 0.04) 
	#q5 = c(0.04, 0.04, 0.04, 0.04, 0.8) 
	#Q = rbind(q1, q2, q3, q4, q5)
	#Q = NULL
	#for (j in 1:npcs) {
	#	q = rep(0.1 / (npcs-1), npcs)
	#	q[j] = 0.9
	#	Q = cbind(Q, q)
	#}
	Q = matrix(0,npcs,q); diag(Q) = 1
	Y = Z %*% Q; N = mvrnorm(n, rep(0, q), diag(q)); Y = Y + N 
	newY = newZ %*% Q; newN = mvrnorm(n, rep(0, q), diag(q)); newY = newY + newN
		
	## Method 2
	#r = svd(X)
	#newr = svd(newX)
	#Y = as.matrix(r$u[ ,1])
	#newY = as.matrix(newr$u[ ,1])
	#for (j in 2:npcs)  {
	#	Y = cbind(Y, r$u[ ,j]) 
	#	newY = cbind(newY, newr$u[ ,j]) 
	#}	
	##Y = Y + mvrnorm(n, rep(0, q), diag(q))
	##newY = newY + mvrnorm(n, rep(0, q), diag(q))

	# Draw heatmap
	#pdf(paste("sim2_hmap.pdf",sep=""))
	#my_palette <- colorRampPalette(c("blue", "yellow"))(n = 299)
	#myHclustDist = function(x) { hclust(d=x, method="average") }
	#heatmap.2(t(X), dendrogram="none", Rowv=F, Colv=F, symm=F, col = my_palette, key=F, keysize=1.5, density.info="none", trace="none", labRow = F)
	#dev.off()
	
	##Measure execution time
	start <- proc.time()
	if (method == "spls") {
		res = spls(X, Y, K=q, eta=0.5, kappa=0.5, select="pls2", fit="simpls", scale.x=TRUE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=FALSE)
	} else if (method == "respls") {
		res = opt(X, Y, L, V0, S, invS, nrmS, q, nlam,outdir)
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
	V_hat2 = matrix(0,p,npcs)	
	V_hat = matrix(0,p,(npcs*nlam))
	for (j in 1:npcs) {
		min_j = 1
		min_ang = 1e6
		for (j2 in 1:ncol(W)) { 
			ang = acos( abs ( crossprod(V0[ ,j], W[ ,j2]) ) / ( norm(V0[ ,j])*norm(W[ ,j2]) ) )
			if (!(ang>=0 & ang<=pi)) 
				stop("[2]Between-vector angle value is wrong!")
			if (ang<min_ang) { 
				min_ang = ang
				min_j = j2				
			}
		}
		V_hat2[, j] =  W[ ,min_j]
			
		#DEBUG
		print(paste("[runParallel3]iter=",i,", j=",j,", min_j=",min_j,", min_ang=",min_ang,sep=""))
		
		if (method == "respls") {
			V_hat[, ((j-1)*nlam + 1):(j*nlam)] =  res$bW[ ,((min_j-1)*nlam + 1):(min_j*nlam)]
		}	
	}

	# Were components correctly split?
	#csplit <- length(unique(best_match)) / length(best_match)
	
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
	v0 <- V0[ ,1]		
	if (npcs>1)
		for (j in 2:ncol(V0)) {
			v0 <- v0 | V0[ ,j]
		}
	#print("v0")
	#print(which(v0!=0))
	w <- W[ ,1]
	if (ncol(W)>1)
		for (j in 2:ncol(W)) {
			w <- w | W[ ,j]
		}
	#print("w")
	#print(which(w!=0))
	t_tp = sum(v0 & w)
	t_fp = sum((!v0) & w)
	t_tn = sum(!v0) - t_fp
	t_prec = t_tp/(t_tp+t_fp)
	
	if (method == "spls") {
		res = list(X=res$X,Y=res$Y,U=U,betapls=res$betapls,elapsed=elapsed,mspe=mspe,recall2=recall2,prec2=prec2,t_prec=t_prec)
	} else if (method == "respls") {
		res <- list(X=res$X,Y=res$Y,U=U,P=res$P,Q=res$Q,betapls=res$betapls,elapsed=elapsed,var=var,cpev=cpev,mspe=mspe,recall=recall,prec=prec,recall2=recall2,prec2=prec2,t_prec=t_prec)
	}
	return(res)

}
