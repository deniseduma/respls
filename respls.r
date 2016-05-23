#source("kernel2.r")
source("proxgrad4.r")

source("mybic.r")
source("respls_v.r")
#source("respls_admm_v.r")
source("respls_admm_w.r")

respls = function(X, Y, L, ncomp, scale.x, scale.y, nlambdas, max.iter, tol, outdir) {
	
	n<-nrow(X)
        p<-ncol(X)
        q<-ncol(Y)

	stop <- FALSE
	number <- ncomp
		
	# Center X
	one <- rep(1, n)
	meanx <- drop( one %*% X ) / n
    	X <- scale( X, meanx, FALSE )
    	# Center Y
	meany <- drop( one %*% Y ) / n
    	Y <- scale( Y, meany, FALSE )

    	# Scale X
	if ( scale.x )
	{
		normx <- sqrt( drop( one %*% (X^2) ) / (n-1) )
        	if ( any( normx < .Machine$double.eps) )
        	{ stop("Some of the columns of the predictor matrix have zero variance.") }
       		 X <- scale( X, FALSE, normx )
	} else { normx <- rep( 1, p ) }
	
	# Scale Y
    	if ( scale.y )
	{
		normy <- sqrt( drop( one %*% (Y^2) ) / (n-1) )
        	if ( any( normy < .Machine$double.eps ) )
        	{ stop("Some of the columns of the response matrix have zero variance.") }
        	Y <- scale( Y, FALSE, normy )
    	} else { normy <- rep( 1, q ) }
	
	rho_v <- 1
	
	YtXs <- info <- list()
	YtX0 <- YtX <- Matrix(crossprod(Y, X))
	# Total variance tO be explained
	tvar <- sum( diag (crossprod(YtX0) ) )

	# Eigendecomp of L used to invert 2 lam2 L + rho I 
	Le <- eigen(L, symmetric=TRUE)

	V <- W <- W2 <- NULL
	for(kk in 1:ncomp){
		
		gc()
		cat("\n\n")
		print(paste("[respls]PC", kk, sep=""))
		
		K <- crossprod(YtX)
		
		# Init v 
		#r <- svd(YtX)
		#v <- Matrix(rW[, 1])
		r <- svd(t(YtX))
		rU <- r$u
		rs <- r$d
		rW <- r$v
		v <- Matrix(rU[, 1])
		print("dim(rU)")
		print(dim(rU))
		print("dim(rW)")
		print(dim(rW))
		print("length(rs)")
		print(length(rs))
		
		idx_lam <- cpev <- 0
		u_v <- u_w <- w <- ws <- us <- temp <- NULL
				
		# Compute MV using the matrix inversion lemma 
		#S <- solve(Diagonal(q) + (1/rho_v)*YtX%*%t(YtX))
		#MV <- (1/rho_v)*Diagonal(p) - (1/rho_v^2)*((t(YtX)%*%S)%*%YtX)
						
		# Matrix inits
		#I_k <- NULL
		#if (!is.null(V)) {
		#	I_k <- Matrix(0, kk-1, kk-1)
		#}
		
		# for v
		#DVM <- init_matrices(MV,V,I_k)
				
		i <- 0
      		diff <- 1
     		while((i<max.iter) & (diff>tol)){
			i<-i+1
			cat("\n")
			print(paste("[respls]v-w iteration ", i,sep=""))
			
			# Update w
			newY <- YtX %*% v
			res <- respls_admm_w(YtX, newY, L, Le, NULL, w, u_w, nlambdas)
			ws <- res$ws
			us <- res$us

			# Normalize columns of ws
			ws <- ws %*% Diagonal(ncol(ws), 1/sqrt(colSums(ws^2)))
			
			# Select the best lambda using BIC 
			res <- select_best_lambda(ws, W, YtX0, YtX, newY, tvar)
			bic <- res$bic
			cpev <- res$cpev
			idx_lam <- res$idx_lam
			
			# Best w selected
			w <- w2 <- as(ws[, idx_lam], "sparseMatrix")
			u_w <- Matrix(us[, idx_lam])
						
			if (nnzero(w)==0) {
				number <- kk-1
				stop <- TRUE
				print(paste("[respls]No non-zero coeffs found for component ",kk,sep=""))
				gc()
				break
			}
									
			# Re-estimate w via OLS
			#sYtX <- Matrix(crossprod(Y, X[, which(w!=0)]))
			sYtX <- Matrix(YtX[, which(w!=0)],q,nnzero(w))
			#sv <-Matrix(v[which(w!=0)],nnzero(w),1)
			#sY <- sYtX %*% sv
			#FIXME
			print("dim(sYtX)")
			print(dim(sYtX))
			print("[respls]cond(crossprod(sYtX))")
			print(kappa(crossprod(sYtX)))
			print(kappa(crossprod(sYtX)+0.001*Diagonal(ncol(sYtX))))
			w_ols <- solve(crossprod(sYtX)+0.001*Diagonal(ncol(sYtX)), crossprod(sYtX,newY))
			w2[which(w!=0)] <- w_ols
			
			# Check stop condition
			normw <- sqrt(sum(w2^2))
			w3 <- w2/normw
			if (!is.null(temp)) {
				diff<-max(abs(w3-temp))
				if (diff>0) {
					aux <- w3-temp
					maxpos=which(abs(aux)==diff)
					print(paste("[respls]maxpos=", maxpos, sep=""))
					print(paste("[respls]max(temp)=",temp[maxpos],", max(w3)=", w3[maxpos],sep=""))
				}
				print(paste("[respls]diff ", diff, sep=""))
			}
			temp <- w3
			
			# Update v
			#newY <- YtX %*% w
			#res <- respls_admm_v(YtX, newY, K, V, NULL, NULL, v, u_v, rho_v)
			res <- respls_v(t(YtX), rU, rs, rW, V, w)
			v <- res$v
			u_v <- res$u
		
		}#end while
		
		if (stop==TRUE) 
			break
				
		if(i==max.iter){
			number <- kk-1
			stop <- TRUE
			cat("[respls]Fail to converge! Increase the number of iterations !","\n")
			gc()
			break
		}
		
		# Update W
		if (is.null(W)) {
			W <- w
		} else {	
			W <- cbind2(W, w)
		}
		
		# Update W2
		#if (is.null(W2)) {
		#	W2 <- w2
		#} else {	
		#	W2 <- cbind2(W2, w2)
		#}
		
		# Update V
		if (is.null(V)) {
			V <- Matrix(v)
		} else {	
			V <-cbind2(V, v)
		}
		
		#DEBUG
		cat("\n")
		print("[respls]VtV")
		print(crossprod(V))
		
		# Return info
		YtXs[[kk]] <- YtX
		info[[kk]] <- list(v=v,w=w,ws=ws,idx_lam=idx_lam,cpev=cpev)

		## Update residual matrix YtX
		YtX <- YtX - YtX%*%(v%*%t(v))
		
		# Remaining var
		Pr <- YtX0%*%(V%*%t(V))
		cvar <- sum(diag(crossprod(Pr))) 
		rvar <- cvar / tvar
		print(paste("tvar ", tvar, ", cvar ", cvar,", rvar ",rvar,sep=""))
		if ((1-rvar)<=0.1)
			break
		
	}#end for k	
	
	# Compute betapls
	U <- X%*%W
	Q <- solve(crossprod(U), crossprod(U,Y))
	betapls <- W%*%Q

	object <- list(ncomp=number,YtXs=YtXs,info=info,betapls=betapls,meanx=meanx,meany=meany,normx=normx,normy=normy)
    	class(object) <- "respls"
    	object
}

#init_matrices <- function(M,Z,I_k) {
#	DZM <- NULL
#	if (!is.null(Z)) { 
#		ZM <- t(Z)%*%M
#		ZMZ <- ZM%*%Z
#		D <- solve(I_k + ZMZ)
#		DZM <- D%*%ZM
#	} 
#	return(DZM)
#}

predict <- function(object, newX) {

	betapls <- object$betapls
	p <- length(betapls)

        if ( ncol(newX)!=p){ stop("The dimension of test dataset is inapproprite!") }
	
        # Predict
        newX <- scale( newX, object$meanx, object$normx )
	pred <- newX %*% betapls + matrix(1,nrow(newX),1) %*% object$meany

	return(pred)
}

