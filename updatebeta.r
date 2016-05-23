#source("sampleTPoints.r")

#Soft thresholding and thourough search for lambda to update beta 

updatebeta <- function(X, K, alpha, p.ini, pcer, thr.range, sampleNo, sampleSize, savepath) {
	
	p.ini  <- p <- ncol(X)
	err <- pcer*p.ini
	stop <- FALSE
	
	ols <- K %*% alpha
	#lambdas <- sort(c(abs(ols),0),decreasing=TRUE)nlam1 <- 30
	nlam1 <- 30
	max_lam1 <- max(abs(ols))
	if (n<p){
		min_lam1 <- 0.01*max_lam1
	} else {
		min_lam1 <- 0.01*max_lam1
	}
	lambdas <- logspace(log10(min_lam1),log10(max_lam1),nlam1)
	lambdas <- sort(c(lambdas,0),decreasing=TRUE)

	cat("\n")
	print("[updatebeta]lambdas")
	print(lambdas[c(1, length(lambdas)-1)])
	
	if(savepath) 
		selprobpath <- matrix(nrow=length(ols),ncol=length(lambdas))
	
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	for(l in 1:length(lambdas)){
		cat("\n")
		print(paste("[updatebeta]l=",l,", lambda=",lambdas[l],sep=""))
		
		temp <- softOp2(X,alpha,lambdas[l]/2,sampleNo,sampleSize)
		t <- temp!=0
		qs[l] <- mean(colSums(t))
		sp <- rowMeans(t)
		
		if(savepath)
			selprobpath[,l] <- sp
		
		thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
		if(thrall[l]>=thr.range[1]){
			stable <- which(sp>=thrall[l])
			if (length(stable)>0) {
				cat("\n")
				print(paste("[updatebeta]Within range, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
				print(paste("[updatebeta]length(stable)=",length(stable),sep=""))
				ls <- l
				break
			}
		}
		
		print(paste("[updatebeta]thr=",thrall[l],", qs=",qs[l],sep=""))
	}#end for
	
	thr <- thrall[ls]
	if(thr > thr.range[2]){
		cat("\n")
		print(paste("[updatebeta]Too big, lambda=",lambdas[ls],", qs=",qs[ls],", thr=",thrall[ls],sep=""))
		while(pcer <= 0.5){
			pcer <- pcer + 0.01	
			thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
			thr <- thrall[ls]
			if(thr < thr.range[2])break
		}
		print(paste("[updatebeta]New pcer=",pcer,", new thr=",thr,sep=""))
	}
	
	stable <- which(sp>=thr)
	print(paste("[updatebeta]length(stable)=",length(stable),sep=""))
	if(length(stable)==0) 
		stop <- TRUE
 	
	lambda = lambdas[ls]/2
	beta <- (sign(ols)*(abs(ols)>=lambda)*(abs(ols)-lambda))
	if(savepath)
		return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=2*lambda,l=ls,selprobpath=selprobpath[, 1:ls]))
	else
		return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=2*lambda,l=ls))
		
}

softOp2 <- function(X, alpha, lambda, sampleNo, sampleSize) {
	#subsets <- sapply(1:sampleNo, function(s) {sampleTPoints(rownames(X))})
	subsets <- sapply(1:sampleNo, function(s) {sample(1:nrow(X), sampleSize*nrow(X))})
	res <- sapply(1:sampleNo,function(s) {soft2(s,subsets,X,alpha,lambda)})
	return(res)
}

soft2 <- function(index,subsets,X,alpha,lambda){
	
	X <- X[subsets[, index], ]
	K <- t(X) %*% X
	
	##rbf <- rbfdot(sigma)
	##K <- kernelMatrix(rbf, t(X))
	
	#X <- matrix(0,nrow(M),ncol(M))	
	#for (i in 1:nrow(M)) 
	#	for (j in 1:ncol(D))
	#		X[i,j] <- rnorm(1,M[i,j],D[i,j])
	#K <- computeKernel(X)
	
	ols <- K %*% alpha
	
	ols <- sign(ols)*(abs(ols)>=lambda)*(abs(ols)-lambda)
	ols[is.na(ols)] <- 0
	return(ols)
}




