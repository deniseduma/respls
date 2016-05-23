library("rdetools")

source("admm_lasso2.r")
source("sampleTPoints.r")

admm2 <- function(X, XtX, Lapl, S, alpha, x_init=NULL, p.ini, nlambdas, sselection, pcer, thr.range, sampleNo, sampleSize, savepath) {
		
	n <- nrow(X)
	p.ini  <- p <- ncol(X)
	if (sselection) {
		err <- pcer*p.ini
		stop <- FALSE
	}
	
	rho <- 1
	relaxpar <- 1.8
	
	x <- Matrix(numeric(p))
	u <- Matrix(numeric(p))
	if(!is.null(x_init)) 
		x <- x_init
	
	max_lam <- max(abs(XtX%*%alpha))
	#if (n>p){
	#	min_lam <- 0.0001*max_lam
	#} else {
	#	min_lam <- 0.01*max_lam
	#}
	lambdas <- logspace(log10(0.01*max_lam),log10(0.95*max_lam),nlambdas)
	lambdas <- sort(abs(lambdas),decreasing=TRUE)
	
	print("[admm]lambdas")
	print(lambdas[c(1, length(lambdas))])
	
	if (!sselection)
		xxs <- Matrix(0, nrow=p, ncol=length(lambdas))
	
	if(savepath) 
		selprobpath <- Matrix(0, nrow=p,ncol=length(lambdas))
	
	if (sselection) {
		qs <- numeric(length(lambdas))
		thrall <- numeric(length(lambdas))
		ls <- length(lambdas)
	}
	
	y <- X %*% alpha
	#y <- Matrix(scale(y, center=F, scale=T))
	
	for(l in 1:length(lambdas)){
		#cat("\n")
		#print(paste('[', l, ']lambda=', lambdas[l],sep=""))
		
		lam1 <- lambdas[l] #0.8*lambdas[l]
		lam2 <- 0.3*lambdas[l]
		
		#cache the factorization
		res <- factor2(XtX,Lapl,lam2,rho)
		L <- res$L 
		U <- res$U
		
		temp <- sss(X,XtX,y,Lapl,L,U,x,u,lam1,lam2,rho,relaxpar,sselection,sampleNo,sampleSize)
		#temp <- (1 + lam2) * temp

		if (!sselection) {
			sp <- temp$x
		} else {
			t <- temp$x!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			stable <- which(sp>=thrall[l])
		}
		
		if (!sselection)
			xxs[,l] <- sp

		if(savepath)
			selprobpath[,l] <- sp
				
		if (sselection) {
			if(thrall[l]>=thr.range[1]){
				#if (length(stable)>0) {
					print(paste("[admm]Within range, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
					print(paste("[admm]length(stable)=",length(stable),sep=""))
					ls <- l
					break
				#}
			}
		}
		
		#warm start
		x <- sp
		u <- temp$u
	}#end for
		
	if (sselection) {
		thr <- thrall[ls]
		if(thr > thr.range[2]){
			print(paste("[admm]Too big, lambda=",lambdas[ls],", qs=",qs[ls],", thr=",thrall[ls],sep=""))
			while(pcer <= 0.5){
				pcer <- pcer + 0.01	
				thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
				thr <- thrall[ls]
				if(thr < thr.range[2])break
			}
			print(paste("[admm]New PCER=", pcer,", new thr=", thr, sep=""))
		}	 
	
		stable <- which(sp>=thr)
		print(paste("[admm]New length(stable)=",length(stable),sep=""))
		if(length(stable)==0) 
			stop <- TRUE
	
		#Solve for beta
		#FIXME solve for beta using the stable vars directly
		x<-lasso(X,y,gamma,lam2,rho,relaxpar,sselection)
	} 
	
	#if (savepath) {
		#plot regularization or selection prob path
		#dev.new()
		#matplot(log10(lambdas), t(selprobpath), cex=1, pch=16, type="l", xaxt="n", xlim=c(log10(min_lam1), log10(max_lam1)), xlab ="lambda", ylab="Coeffs/Selection prob", main="Regularization path", cex.lab=1.25, cex.axis=1.25)
		#Sys.sleep(0.5) 
		#dev.off()
	#}
		
	#if (sselection) {
	#	if(savepath)
	#		return(list(beta=beta,sp=as.vector(sp),thr=thr,stop=stop,lambda=lambdas[ls],l=ls,selprobpath=selprobpath[, 1:ls]))
	#	else
	#		return(list(beta=beta,sp=as.vector(sp),thr=thr,stop=stop,lambda=lambdas[ls],l=ls))
	#} else {	
		if (savepath)
			return(list(beta=xxs,thr=(-Inf),selprobpath=selprobpath))
		else
			return(list(beta=xxs,thr=(-Inf)))
	#}	

}	

sss <- function(A,AtA,b,Lapl,L,U,x,u,lam1,lam2,rho,relaxpar,sselection,sampleNo,sampleSize) {
	if (sselection) {
		#subsets <- sapply(1:sampleNo, function(s) {sampleTPoints(rownames(X))})
		#res <- sapply(1:sampleNo, function(index) {lasso(X[subsets[, index],],y,lambda,rho,relaxpar,sselection)})
	} else {
		res <- admm_lasso2(A,AtA,b,Lapl,L,U,x,u,lam1,lam2,rho,relaxpar,sselection)
	}
	return(res)
}


