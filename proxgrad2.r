library("rdetools")
source("sampleTPoints.r")

proxgrad <- function(X, K, L, alpha, init_beta=NULL, p.ini, tol, sselection, pcer, thr.range, sampleNo, sampleSize, savepath) {

	n <- nrow(X)
	p.ini  <- p <- ncol(X)
	if (sselection) {
		err <- pcer*p.ini
		stop <- FALSE
	}
		
	nlam1 <- 30
	max_lam1 <- max(abs(K%*%alpha))
	max_lam1 <- max_lam1/0.8
	if (n<p){
		min_lam1 <- 0.05*max_lam1
	} else {
		min_lam1 <- 0.01*max_lam1
	}
	lambdas <- logspace(log10(min_lam1),log10(max_lam1),nlam1)
	lambdas <- sort(abs(lambdas),decreasing=TRUE)
	
	print("[proxgrad]lambdas")
	print(lambdas[c(1, length(lambdas))])
	#print(lambdas)

	if(is.null(init_beta))
		init_beta <- numeric(p.ini)

	if(savepath) 
		selprobpath <- matrix(nrow=p,ncol=length(lambdas))
	
	if (sselection) {
		qs <- numeric(length(lambdas))
		thrall <- numeric(length(lambdas))
		ls <- length(lambdas)
	}
	
	for(l in 1:length(lambdas)){
		cat("\n")
		print(paste('[', l, ']lambda=', lambdas[l],sep=""))
		
		temp <- sss(X,L,alpha,init_beta,0.8*lambdas[l],0.2*lambdas[l],tol,sselection,sampleNo,sampleSize)
		if (sselection==F) {
			sp <- temp 
		} else {
			t <- temp!=0
			qs[l] <- mean(colSums(t))
			sp <- rowMeans(t)
			thrall[l] <- ((qs[l]^2/(err*p.ini))+1)/2
			stable <- which(sp>=thrall[l])
		}
		
		if(savepath)
			selprobpath[,l] <- sp

		if (sselection) {
			if(thrall[l]>=thr.range[1]){
				#if (length(stable)>0) {
					print(paste("[proxgrad]Within range, lambda=",lambdas[l],", qs=",qs[l],", thr=",thrall[l],sep=""))
					print(paste("[proxgrad]length(stable)=",length(stable),sep=""))
					ls <- l
					break
				#}
			}
		}
		
		#warm start
		init_beta <- sp
	}#end for

	if (sselection) {
		thr <- thrall[ls]
		if(thr > thr.range[2]){
			print(paste("[proxgrad]Too big, lambda=",lambdas[ls],", qs=",qs[ls],", thr=",thrall[ls],sep=""))
			while(pcer <= 0.5){
				pcer <- pcer + 0.01	
				thrall <- ((qs^2/((pcer*p.ini)*p.ini))+1)/2
				thr <- thrall[ls]
				if(thr < thr.range[2])break
			}
			print(paste("[proxgrad]New PCER=", pcer,", new thr=", thr, sep=""))
		}	 
	
		stable <- which(sp>=thr)
		print(paste("[proxgrad]New length(stable)=",length(stable),sep=""))
		if(length(stable)==0) 
			stop <- TRUE
	
		#Solve for beta
		#FIXME solve for beta using the stable vars directly
		beta<-regression_graph(X,L,alpha,init_beta,0.8*lambdas[ls],0.2*lambdas[ls],tol,sselection)
	} else {
		#FIXME; if no stability selection
		t <- selprobpath!=0
		#print("[proxgrad]selprobpath")
		#print(selprobpath)
		beta <- rowMeans(t)
		#print("[proxgrad]beta")
		#print(t(beta))
		beta[beta<thr.range[1]] <- 0
		print(paste("[proxgrad]nnzero(beta)=", (nnzero(beta)/p)*100,sep=""))
		#X <- X[, which(beta>0)]
		#y <- X %*% alpha[which(beta>0)]
		#small <- solve((t(X)%*%X))%*%t(X)%*%y
		#beta[beta>0] <- small
	}
	
	#if (savepath) {
		#plot regularization or selection prob path
		#dev.new()
		#matplot(log10(lambdas), t(selprobpath), cex=1, pch=16, type="l", xaxt="n", xlim=c(log10(min_lam1), log10(max_lam1)), xlab ="lambda", ylab="Coeffs/Selection prob", main="Regularization path", cex.lab=1.25, cex.axis=1.25)
		##Sys.sleep(0.5) 
		##dev.off()
	#}
		
	if (sselection) {
		if(savepath)
			return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=lambdas[ls],l=ls,selprobpath=selprobpath[, 1:ls]))
		else
			return(list(beta=beta,sp=sp,thr=thr,stop=stop,lambda=lambdas[ls],l=ls))
	} else {	
		if (savepath)
			return(list(beta=beta,thr=(-Inf),selprobpath=selprobpath))
		else
			return(list(beta=beta,thr=(-Inf)))
	}	

}	

sss <- function(X,L,alpha,init_beta,lam1,lam2,tol,sselection,sampleNo,sampleSize) {
	if (sselection) {
		subsets <- sapply(1:sampleNo, function(s) {sampleTPoints(rownames(X))})
		res <- sapply(1:sampleNo, function(index) {regression_graph(X[subsets[, index],],L,alpha,init_beta,lam1,lam2,tol,sselection)})
	} else {
		res <- regression_graph(X,L,alpha,init_beta,lam1,lam2,tol,sselection)
	}
	return(res)
}

regression_graph <- function(X,L,alpha,init_beta,lam1,lam2,tol,sselection) {
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	y <- X %*% alpha
	beta <- init_beta

	#fixed step size
	#r <- eigen(t(X) %*% X + 2*lam2*L, only.values = TRUE)
	#lip <- max(r$values)
	#t <- 1/lip
	
	if (sselection==F) {
		k2 <- 1
		obj <- numeric()
		obj[k2] <- reg_obj(y,X,beta,L,lam1,lam2)
	}
	
	t<-1
	b<-0.5
	diff <- 1
	crt_obj <- smth_obj(y,X,beta,L,lam1,lam2)
	while (diff>tol) {
		#line search to find step size
		k1 <- 1
		#cat("\n")
		while (TRUE) {
			g <- gprox(y,X,beta,L,lam1,lam2,t)
			new_beta <- beta - t * g
			new_obj <- smth_obj(y,X,new_beta,L,lam1,lam2)
			rhs <- crt_obj - t * t(smth_grad(y,X,beta,L,lam1,lam2)) %*% g + (t/2) * sum(g^2)
			#print(paste("diff=", sqrt(sum(g^2)), ", obj=",new_obj,sep=""))
			if (new_obj <= rhs)	{
				beta <- new_beta
				crt_obj <- new_obj
				break
			}
			t <- b * t
			k1 <- k1 + 1
		}#end while
		#print(paste("[regression_graph]k1=", k1, ", t=", t, sep=""))
	
		diff <- sqrt(sum(g^2))
		#print(paste("[regression_graph]diff=", diff, sep=""))
		
		if (sselection==F) {
			k2 <- k2 + 1
			obj[k2] <- reg_obj(y,X,beta,L,lam1,lam2)
		}
	
	}#end while
	print(paste("[regression_graph]k2=", k2, ", t=", t, sep=""))
	
	#if (sselection==F) {
	#	dev.new()
	#	plot(x<-1:k2, y<-log10(obj), type="p", cex=1, col="red", xlab="iteration", ylab="objective", xlim=c(1,k2), main=paste("lambda ", lam1, sep=""))
#		#Sys.sleep(0.5) 
#		#dev.off()
	#}
	
	return (beta)
}


gprox <- function(y,X,beta,L,lam1,lam2,t) {
	smth_step <- beta - t*smth_grad(y,X,beta,L,lam1,lam2)
	ret <- (1/t) * (beta - proxl1(smth_step,t,lam1))
	return(ret)
}

smth_grad <- function(y,X,beta,L,lam1,lam2) {
	ret <- -t(X) %*% y + t(X) %*% X %*% beta + 2*lam2*L%*%beta
	return(ret)
}

reg_obj <- function(y,X,beta,L,lam1,lam2) {
	obj <- (1/2)*sum((y - X%*%beta)^2) + lam1*sum(abs(beta)) + lam2*t(beta)%*%L%*%beta
	return(obj)
}

smth_obj <- function(y,X,beta,L,lam1,lam2) {
	obj <- (1/2)*sum((y - X%*%beta)^2) + lam2*t(beta)%*%L%*%beta
	return(obj)
}

proxl1 <- function(beta, t, lam) {
	scle <- lam*t
	beta[beta<scle & beta>(-scle)] <- 0
	beta[beta>scle] <- beta[beta>scle] - scle  
	beta[beta<(-scle)] <- beta[beta<(-scle)] + scle  
	return (beta)
}

