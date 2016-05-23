library(MASS)
library(Matrix)
library(rdetools)
library(glmnet)
source("update_v.r")
source("update_w.r")

DEBUG = 0

#rad2deg <- function(rad) {(rad * 180) / (pi)}
#deg2rad <- function(deg) {(deg * pi) / (180)}
angle = (45 * pi) / (180)
maxit = 100
eps = 1e-3
rho_v = 1
rho_w = 1

alpha1 = 0.6
alpha2 = 1-alpha1

opt2 = function(X,Y,L,V0,S,invS,nrmS,npcs,nlam,outdir) {

n = nrow(X)
p = ncol(X)
q = ncol(Y)

# Center X
one = rep(1, n)
meanx = drop(one %*% X) / n
names(meanx) = NULL
X = scale( X, meanx, F )
# Center Y
meany = drop(one %*% Y) / n
names(meany) = NULL
Y = scale( Y, meany, F )
# Scale X
normx = sqrt( drop( one %*% (X^2) ) / (n-1) )
names(normx) = NULL
X = scale(X, F, normx)
# Scale Y
normy = sqrt( drop( one %*% (Y^2) ) / (n-1) )
names(normy) = NULL
Y = scale( Y, F, normy )
XO=X; YO=Y

## New data matrix
X = X0 = (t(Y)%*%X) / (n-1)
XtX = t(X)%*%X
## Total variance to be explained
tvar = sum( diag ( XtX ) )
## Variance by PCs
r = svd(X);d2 = r$d^2;sd2 = sum(d2);ld2 = length(d2)
var = numeric(ld2)
for (j in 1:ld2)
	var[j] = sum(d2[1:j]) / sd2
	
V = W = NULL
bW = matrix(0, p, npcs*nlam)
cpev = numeric(npcs)
	
## Speedup w update
#M_w = L + diag(p);res = factor_w(M_w);L_w0 = res$L ;U_w0 = res$U
#M_w = as.matrix(nearPD(L)$mat);res = factor_w(M_w);L_w0 = res$L ;U_w0 = res$U
#M_w = L;res = factor_w(M_w);L_w0 = res$L;U_w0 = res$U
	
for (k in 1:npcs) {
	cat("\n\n"); print(paste("PC",k,sep=""))
	
	## Init v from SVD of data matrix
	r = svd(X);v = r$v[, 1]
	y = X%*%v; Xty = t(X)%*%y
	print("v1")
	print(v[which(V0[,1]!=0)])
	print("v2")
	print(v[which(V0[,2]!=0)])
	
	## Search for max_lam
	max_lam = max(abs(Xty))
	Cor = Xty / max(1, q-1) 
	indx = which(abs(Cor)==max(abs(Cor))) 
	print(paste("max_lam1 ",max_lam,", indx ",indx,sep=""))
	for( i in seq(Cor[indx],0,length.out=50) ) {
		print(paste("i ",i,sep=""))
		tmp = numeric(p); tmp[indx]=i
		lls=abs(Xty - XtX%*%tmp - alpha2*( max_lam/alpha1 )*L%*%tmp)
		if( lls[indx] > max(lls[-indx]) ) {
			max_lam = lls[indx] 
			print("t(lls)"); print(t(lls))
			#print(paste("max_l ", max_l, sep=""))
			break
		}
	}
	print(paste("max_lam2 ",max_lam,sep=""))
	#max_lam = max_lam / alpha1
	#max_lam = max(abs(Xty)) / alpha1
	## Generate lambda sequence
	lambdas = exp(seq(from=log(max_lam),to=log(0.01*max_lam),length.out=nlam))
	print("lambdas")
	print(lambdas[c(1, length(lambdas))])

	# Speedup w update
	#Z = t(YtX)
	#my_MZ = backsolve(U_w0, forwardsolve(L_w0, Z))
						
	# Speedup v update
	## Method 1: Cholesky decomposition
	#M_v = K + rho_v*diag(p)
	#res = factor_v(M)
	#L_v = res$L 
	#U_v = res$U
	## Method 2: Woodbury formula
	#Vhat <- D <- NULL
	#invM_v = diag(p) - t(YtX) %*% solve( diag(q) + YtX %*% t(YtX) ) %*% YtX
	#if (!is.null(V)) {
	#	Vhat = invM_v %*% V
	#	D = t(Vhat) %*% V
	#}
	## Method 3: Woodbury formula better computation 
	MV = VMV = NULL
	if (!is.null(V)) {
		W_v = solve( diag(q) + X%*%t(X), X%*%V )
		MV = V - t(X)%*%W_v
		VMV = t(MV)%*%V
	}
	
	iter=1; #indvw = 1
	w = u_w = u_v = rep(0, p); nnz_w = 0; err = dfterm = 1e6
	while (iter<maxit) {
		oldw = w; w = rep(0, p)
		bicw = rep(1e6, nlam)
		ws = matrix(0, nrow=p, ncol=nlam)
		cat("\n"); print(paste("iter ",iter,", BIC",sep=""))
		for (j in 1:nlam) { 
			onnz_w = nnz_w
			lam = lambdas[j];lam1 = alpha1*lam;lam2 = alpha2*lam  
			#L_w = sqrt(lam2) * L_w0; U_w = sqrt(lam2) * U_w0
			## FIXME
			U_w = chol(lam2*L + rho_w*diag(p)); L_w = t(U_w)

			## Speedup w update
			#MZ = (1/lam2) * MZ
			#ZMZ = diag(q) +  t(MZ)%*%Z

			## Compute w using ADMM
			#if (is.null(oldw)) {wi = w} else {wi = oldw}
    			res=update_w(X, y, Xty, XtX, L, S, invS, nrmS, L_w, U_w, NULL, NULL, w, NULL, rho_w, lam1, lam2) 
			x_w = res$x; w2 = w = res$z; u_w = res$u; ws[ ,j] = w
			
			#res2 = glmnet(X,Y,alpha=0.999,intercept=F,standardize=F,family="gaussian",nlambda=50)
			#print(colSums(res2$beta!=0))
			if (norm(w)==0)
				next

			## Compute BIC for w
			idx_w = which(w!=0); nnz_w=length(idx_w)
			sX = X[ ,idx_w,drop=FALSE]
			#if (nnz_w!=onnz_w) { 
				## OLS re-estimation of w
				sw = solve ( crossprod(sX) + 0.001*diag(nnz_w), t(sX)%*%y); w2[idx_w] = sw
				err = (norm(y - X%*%w2)^2) / (q)
				#err = (norm(y - X%*%x_w)^2) / (q)
			#}	
			## df estimation
			## Direct matrix inversion
			DF = sX %*% solve( crossprod(sX) + lam2*L[idx_w, idx_w] + lam2*diag(nnz_w), t(sX))
			#DF = sX %*% solve( crossprod(sX) + lam2*L[idx_w, idx_w], t(sX))
			## Cholesky decomposition
			#M_lam = crossprod(sYtX) + lam2*L[idx_w, idx_w] + lam2*diag(nnz_w)
			#res = factor_w(M_lam); L_lam = res$L; U_lam = res$U
			#DF = backsolve(U_lam, forwardsolve(L_lam, t(sYtX)))
			#DF = sYtX %*% DF
			df =  sum(diag(DF)); dfterm = (log(q)/max(1, q-1))*df
			#plty = 2*0.5*log(choose(p, nnz_w))
			bicw[j] = log(err) + dfterm #+ plty
			#if (DEBUG)  
				print(paste("[opt]j=",j,", nnz_w=",nnz_w,", err=",err,", log(err)=",log(err),", df=",dfterm, ", bic=",bicw[j],sep=""))
			

		}# end for j
		
		## Find w corresponding to min BIC
		idx_bic = which.min(bicw)
		w = ws[ ,idx_bic]; #nnz_w = sum(w!=0)
		#if (DEBUG)
			print(paste("best j=",idx_bic,", best nnz=",sum(w!=0),", best bic=",bicw[idx_bic], sep=""))
		
		## Test stop condition on w 
		#if (!is.null(oldw)) { 
			indw = norm(w - oldw) / norm(oldw)
			#if (DEBUG)
				cat("\n"); print(paste("[opt]indw ",indw,", nnz(w) ",sum(w!=0),", norm(w) ",norm(w),sep=""))
			if (indw<eps) 
				break
		#}
				
		## Update v using ADMM
		y = X%*%w
		Xty = t(X)%*%y
		res=update_v(X, y, Xty, XtX, V, MV, VMV, v, NULL, rho_v)
		x_v = res$x; v = res$z; u_v = res$u
				
		## Update y
		y = X%*%v
		Xty = t(X)%*%y
		
		#indvw = norm(v-oldv)/norm(oldv) + norm(w-oldw)/norm(oldw)
		iter = iter + 1
	} #end while
	
	if (DEBUG) { 
		cat("\n")
		print(paste("[opt]k=",k,", iter=",iter,", indw ",indw,", nnz(w)=",sum(which(w!=0)),sep=""))
		print(paste("[opt]norm(w)=",norm(w),", norm(v)=", norm(v),", norm(v-w)=",norm(v-w),sep=""))
	}

	## Update V
	V = cbind(V, v)
	## Update W  
	w = w / norm(w)
	W = cbind(W, w)
	## Update bW
	ws = ws %*% diag(sqrt(colSums(ws^2)))
	bW[ ,((k-1)*nlam + 1):(k*nlam)] = ws
	## Compute cpev
	pW = W %*% ginv( crossprod(W) ) %*% t(W) 
 	cpev[k] = sum( diag( crossprod(X0%*%pW) ) ) / tvar
		
	## Update data matrix
	prV = V %*% solve( crossprod(V) ) %*% t(V) 
	X = X0 - X0 %*% prV
	XtX = crossprod(X)
	
} #end for k
	
# Compute betapls
TT = XO %*% W
P = ginv(crossprod(TT)) %*% t(TT) %*% XO 
Q = ginv(crossprod(TT)) %*% t(TT) %*% YO
betapls = W%*%Q
	
# DEBUG
cat("\n"); print("VtV")
print(round(crossprod(V), 4))
	
## Plot reg paths
for (k in 1:ncol(W)) { 
	
	ws = bW[, ((k-1)*nlam + 1):(k*nlam)]	
	coefs = t(ws)
	coefs = scale(coefs, FALSE, 1/normx)	
	low = min(coefs) - 0.01 * abs(min(coefs))
	up = max(coefs) + 0.01 * abs(max(coefs))
	ylimit = c(low, up)
	
	#dev.new()
	pdf(paste(outdir, "pc",k,"_reg_path.pdf",sep=""))
	plot(1:nlam, 1:nlam, xlab="Steps", ylab="Standardized coefficients", ylim=ylimit,  type="n", main = paste("Regularization path PC",k,sep=""))
	for (l in 1:ncol(coefs)) {
		if (sum(coefs[,l]!=0)!=0)
			lines(1:nlam, coefs[,l], col = l, lty = 1)
	}
	idx = which(coefs[nlam,]!=0)
	axis(4, at = coefs[nlam, idx], labels = colnames(X0)[idx], cex = 0.1, adj=0, las=3, mar=c(8,8,1,4))
	abline(h = 0, lty = 3)
	dev.off()

	# write non-zero selected coeffs to file
	idx = which(W[ ,k]!=0)
	writeLines(colnames(X0)[idx], paste(outdir, "pc",k,"_coeffs.txt",sep=""))
}

object <- list(X=XO,Y=Y0,V=V,W=W,bW=bW,P=P,Q=Q,betapls=betapls,var=var,cpev=cpev,meanx=meanx,normx=normx,meany=meany,normy=normy)
class(object) <- "respls"
return(object)

}

predict.respls <- function(obj, newX) {
	
	#print("[predict.respls]dim(newX)")
	#print(dim(newX))

	X <- obj$X
	betapls <- obj$betapls
	p <- ncol(X)
	
        if (is.null(newX)) {
		pred <- X %*% betapls #+ matrix(1,nrow(X),1) %*% obj$meany
		pred <- scale(pred, F, 1/obj$normy)
		pred <- scale(pred, -obj$meany, F)
	} else {	
        	if ( ncol(newX)!=p){ stop("The dimension of test dataset is inapproprite!") }
		newX <- scale( newX, obj$meanx, obj$normx )
		pred <- newX %*% betapls #+ matrix(1,nrow(newX),1) %*% object$meany
		pred <- scale(pred, F, 1/obj$normy)
		pred <- scale(pred, -obj$meany, F)
	}
	
	return(pred)
}

