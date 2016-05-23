#ADMM splitting by features with lasso trick

source("admm_lasso_multi2.r")
source("sampleTPoints.r")

NCORES <- 1 
SSIZE <-0

admm <- function(X, K, Lapl, S, alpha, init_beta=NULL, p.ini, sselection, pcer, thr.range, sampleNo, sampleSize, savepath) {
		
	n <- nrow(X)
	p.ini  <- p <- ncol(X)
	
	if (sselection) {
		err <- pcer*p.ini
		stop <- FALSE
	}
	
	while (p%%NCORES!=0) {
		NCORES <<- NCORES-1
	}	
	SSIZE <<- p/NCORES
	print(paste("[lasso]NCORES ", NCORES, sep=""))
	
	rho <- 1
	relaxpar <- 1.8

	beta <- init_beta
	if (is.null(beta))
		beta <- Matrix(0, nrow=SSIZE, ncol=NCORES)

	nlam1 <- 100
	max_lam1 <- max(abs(K%*%alpha))
	#max_lam1 <- max_lam1/0.8
	if (n>p){
		min_lam1 <- 0.0001*max_lam1
	} else {
		min_lam1 <- 0.001*max_lam1
	}
	lambdas <- logspace(log10(min_lam1),log10(max_lam1),nlam1)
	lambdas <- sort(abs(lambdas),decreasing=TRUE)
	print("[admm]lambdas")
	print(lambdas[c(1, length(lambdas))])
	
		
	if (!sselection)
		betas <- Matrix(0, nrow=p, ncol=length(lambdas))
	
	if(savepath) 
		selprobpath <- Matrix(0, nrow=p,ncol=length(lambdas))
	
	if (sselection) {
		qs <- numeric(length(lambdas))
		thrall <- numeric(length(lambdas))
		ls <- length(lambdas)
	}
	
	X <- Matrix(X)
	S <- Matrix(S)
	alpha <- Matrix(alpha)
	
	#extended y
	y <- X %*% alpha
	y <- rbind2(y, Matrix(numeric(p)))

	for(l in 1:length(lambdas)){
		cat("\n")
		print(paste('[', l, ']lambda=', lambdas[l],sep=""))
		
		lam1 <- 0.8*lambdas[l]
		lam2 <- 0.2*lambdas[l]
		#lam1 <- lam1/(n+p)
		#lam2 <- lam2/(n+p)

		#extended X
		newX <- (1/sqrt(1+lam2)) * rbind2(X, sqrt(lam2)*t(S))
		
		#modified l1 penalty
		lam1 <- lam1 / sqrt(1+lam2)
	
		temp <- sss(newX,y,beta,lam1,lam2,rho,relaxpar,sselection,sampleNo,sampleSize)
		
		temp <- temp / sqrt(1+lam2)
		
		if (sselection==F) {
			sp <- temp
		} else {
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			stable <- which(sp>=thrall[l])
		}
		
		if (!sselection)
			betas[,l] <- sp

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
		beta <- sp
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
		beta<-lasso(X,y,gamma,lam2,rho,relaxpar,sselection)
	} #else {
		#FIXME; if no stability selection
		#t <- (selprobpath!=0)
		#beta <- rowMeans(t)
		#beta[beta<thr.range[1]] <- 0
		#print(paste("[admm]nnzero(beta)=", nnzero(beta), sep=""))
		##X <- X[, which(beta>0)]
		##y <- X %*% alpha[which(beta>0)]
		##small <- solve((t(X)%*%X))%*%t(X)%*%y
		##beta[beta>0] <- small
	#}
	
	#if (savepath) {
		#plot regularization or selection prob path
	#	dev.new()
	#	matplot(log10(lambdas), t(selprobpath), cex=1, pch=16, type="l", xaxt="n", xlim=c(log10(min_lam1), log10(max_lam1)), xlab ="lambda", ylab="Coeffs/Selection prob", main="Regularization path", cex.lab=1.25, cex.axis=1.25)
		#Sys.sleep(0.5) 
		#dev.off()
	#}
		
	if (sselection) {
		if(savepath)
			return(list(beta=beta,sp=as.vector(sp),thr=thr,stop=stop,lambda=lambdas[ls],l=ls,selprobpath=selprobpath[, 1:ls]))
		else
			return(list(beta=beta,sp=as.vector(sp),thr=thr,stop=stop,lambda=lambdas[ls],l=ls))
	} else {	
		if (savepath)
			return(list(beta=betas,thr=(-Inf),selprobpath=selprobpath))
		else
			return(list(beta=betas,thr=(-Inf)))
	}	

}	

sss <- function(A,b,xold,lam1,lam2,rho,relaxpar,sselection,sampleNo,sampleSize) {
	if (sselection) {
		#subsets <- sapply(1:sampleNo, function(s) {sampleTPoints(rownames(X))})
		#res <- sapply(1:sampleNo, function(index) {lasso(X[subsets[, index],],y,lambda,rho,relaxpar,sselection)})
	} else {
		res <- admm_lasso_multi2(A,b,xold,lam1,rho,relaxpar,sselection,NCORES,SSIZE)
	}
	return(res)
}


