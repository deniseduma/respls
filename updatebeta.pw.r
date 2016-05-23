library("rdetools")
source("sampleTPoints.r")

##MS
updatebeta.pw <- function(X, K, alpha, p.ini, pcer, thr.range, sampleNo, sampleSize, l=NULL) {
	
	n <- nrow(X)
	p.ini  <- p <- ncol(X)
	err <- pcer*p.ini
	stop <- FALSE
	
	##
	#K <- t(X) %*% X
	#rbf <- rbfdot(sigma)
	#K <- kernelMatrix(rbf, t(X))
	
	ols <- K %*% alpha
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	#nlam1 <- 30
	#max_lam1 <- max(abs(ols))
	#if (n<p){
	#	min_lam1 <- 0.01*max_lam1
	#} else {
	#	min_lam1 <- 0.01*max_lam1
	#}
	#lambdas <- logspace(log10(min_lam1),log10(max_lam1),nlam1)
	#lambdas <- sort(c(lambdas,0),decreasing=TRUE)
	cat("\n")
	print("[updatebeta.pw]lambdas")
	print(lambdas[c(1, length(lambdas)-1)])
	
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	
	if(is.null(l)) l <- which(lambdas==quantile(lambdas,0.5,type=1))[1]
	
	#search for a lambda
	l.min <- 1
	l.max <- length(lambdas)
	for(g in 1:length(lambdas)){
		cat("\n")
		print(paste("[updatebeta.pw]l=",l,", lambda=",lambdas[l],sep=""))
		
		temp <- softOp(X,alpha,lambdas[l]/2,sampleNo,sampleSize)
		
		t <- temp!=0
		qs[l] <- mean(colSums(t))
		sp <- rowMeans(t)
		thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
		stable <- which(sp>=thrall[l])
		
		if(thrall[l] >= thr.range[1] & thrall[l] <= thr.range[2]){
			print(paste("[updatebeta.pw]Within range, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
			print(paste("[updatebeta.pw]length(stable)=",length(stable),", max(sp)=",max(sp),sep=""))
			#if (length(stable)>0) {
				ls <- l
				break
			#}	
		}
		
		if(thrall[l] < thr.range[1]){
			print(paste("[updatebeta.pw]Too small, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
			print(paste("[updatebeta.pw]length(stable)=",length(stable),", max(sp)=",max(sp),sep=""))
			
			l.min <- l
			if(l == length(lambdas))break
			if(thrall[l+1]> thr.range[2]){
				ls <- l+1
				if(thrall[l+1]>1) ls <-l
				
				cat("\n")
				print(paste("[updatebeta.pw]ls=",ls,sep=""))
				temp <- softOp(X,alpha,lambdas[ls]/2,sampleNo,sampleSize)
				t <- temp!=0
				qs[ls] <- mean(colSums(t))
				sp <- rowMeans(t)
				thrall[ls] <- ((qs[ls]^2/(err*p.ini))+1)/2
				break
			}
			l <- min(length(lambdas),l.max,l + ceiling(length(lambdas)/(g+1)))
			while(thrall[l]!=0 ){  # if thr for current lambda available decrease lambda 
				l <- l-1
				if(l == 0)break
			} 
		}
			
		if(thrall[l] > thr.range[2]){ 
			print(paste("[updatebeta.pw]Too big, lambda=", lambdas[l], ", qs=",qs[l],", thr=",thrall[l],sep=""))
			print(paste("[updatebeta.pw]length(stable)=",length(stable),", max(sp)=",max(sp),sep=""))
			
			l.max <- l
			if(l == 0)break
			if(thrall[l-1]<thr.range[1]&thrall[l-1]!=0){
				ls <- l
				break
			}
			l <- max(1,l.min,l - ceiling(length(lambdas)/(g+1)))
			while(thrall[l]!=0){ 
				l <- l+1
				if(l == length(lambdas))break
			} 
		}
	}#end for
	
	thr <- ((qs[ls]^2/((pcer*p.ini)*p.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0) 
		stop <- TRUE
	
	lambda = lambdas[ls]/2
	beta <- (sign(ols)*(abs(ols)>=lambda)*(abs(ols)-lambda))
	return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=2*lambda,l=ls))
}

##DD
updatebeta.pw2 <- function(X, alpha, pcer, p.ini, thr.range, sampleNo, sampleSize,l=NULL) {
	p.ini  <- p <- ncol(X)
	err <- pcer*p.ini
	stop <- FALSE
	
	##
	K <- t(X) %*% X
	ols <- K %*% alpha
	##
	#with rbf
	#ols <- X %*% alpha
	
	lambdas <- sort(c(abs(ols),0),decreasing=TRUE)
	cat("\n")
	print("[updatebeta.pw]lambdas")
	print(lambdas[c(1, length(lambdas)-1)])
	
	qs <- numeric(length(lambdas))
	thrall <- numeric(length(lambdas))
	ls <- length(lambdas)
	
	if(!is.null(l)) {
		cat("\n")
		print(paste("[updatebeta.pw]l=",l,sep=""))
		
		temp <- softOp(X,alpha,lambdas[l]/2,sampleNo,sampleSize)
		t <- temp!=0
		qs[l] <- mean(colSums(t))
		sp <- rowMeans(t)
		thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
		
		if(thrall[l] >= thr.range[1] & thrall[l] <= thr.range[2]){
			print(paste("within range, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
			stable <- which(sp>=thrall[l])
			print(paste("[updatebeta.pw]length(stable)=",length(stable),sep=""))
			
			ls <- l
			
			thr = thrall[ls]
			stable <- which(sp>=thr)
			if(length(stable)>0) {
				lambda = lambdas[ls]/2
				beta <- (sign(ols)*(abs(ols)>=lambda)*(abs(ols)-lambda))
				return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=2*lambda,l=ls))
			}	
		}
	}
	
	#binary search for a lambda
	l.min <- 1
	l.max <- length(lambdas)
	while (l.max >= l.min) {
		#l <- which(lambdas==quantile(lambdas[l.min:l.max],0.5,type=1))[1]
		l <- ceiling((l.max - l.min + 1)/2)

		cat("\n")
		print(paste("[updatebeta.pw]l=",l,sep=""))
				
		temp <- softOp(X,alpha,lambdas[l]/2,sampleNo,sampleSize)
		t <- temp!=0
		qs[l] <- mean(colSums(t))
		sp <- rowMeans(t)
		thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
		
		if(thrall[l] >= thr.range[1] & thrall[l] <= thr.range[2]){
			print(paste("within range, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
			stable <- which(sp>=thrall[l])
			print(paste("[updatebeta.pw]length(stable)=",length(stable),sep=""))
			ls <- l
			break
		}
		
		if(thrall[l] < thr.range[1]) {
			print(paste("too small, lambda=",lambdas[l], ", qs=",qs[l],", thr=",thrall[l],sep=""))
			stable <- which(sp>=thrall[l])
			print(paste("[updatebeta.pw]length(stable)=",length(stable),sep=""))
			
			l.min <- l + 1 
		}	

		if(thrall[l] > thr.range[2]) { 
			print(paste("too big, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
			stable <- which(sp>=thrall[l])
			print(paste("[updatebeta.pw]length(stable)=",length(stable),sep=""))
			
			l.max <- l - 1
		}	
		

	}#end while	
	
	if (l.max<l.min) 
		print(paste("[updatebeta.pw]No suitable lambda found, ls=", ls,", lambda=",lambdas[ls],sep=""))
	
	thr <- ((qs[ls]^2/((pcer*p.ini)*p.ini))+1)/2
	stable <- which(sp>=thr)
	if(length(stable)==0) 
		stop <- TRUE
	
	lambda = lambdas[ls]/2
	beta <- (sign(ols)*(abs(ols)>=lambda)*(abs(ols)-lambda))
	return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=2*lambda,l=ls))
}

softOp <- function(X, alpha, lambda, sampleNo, sampleSize) {
	#subsets <- sapply(1:sampleNo, function(s) {sampleTPoints(rownames(X))})
	subsets <- sapply(1:sampleNo, function(s) {sample(1:nrow(X), sampleSize*nrow(X))})
	res <- sapply(1:sampleNo,function(s) {soft(s,subsets,X,alpha,lambda)})
	#res <- sapply(1:sampleNo,function(s) {soft(s,M,D,alpha,lambda)})
	return(res)
}

soft <- function(index,subsets,X,alpha,lambda){
	
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

soft0<-function(a,para){
  b<-sort(abs(a))
  b<-abs(a)-para
  b<-(b+abs(b))/2
  b<-sign(a)*b
  b
}


