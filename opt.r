library(MASS)
library(pls)
library(rdetools)
source("update_v.r")
source("update_w.r")

DEBUG = 0

#rad2deg <- function(rad) {(rad * 180) / (pi)}
#deg2rad <- function(deg) {(deg * pi) / (180)}
angle = (45 * pi) / (180)
maxit = 100
eps = 1e-2
rho_v = 1
rho_w = 1

#FIXME
alpha1 = 0.2 #0.6
alpha2 = 1 - alpha1

opt = function(X,Y,L,V0,S,iS,nS,npcs,nlam,outdir) {

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
X = scale( X, F, normx )
# Scale Y
normy = sqrt( drop( one %*% (Y^2) ) / (n-1) )
names(normy) = NULL
Y = scale( Y, F, normy )

##FIXME
## New data matrix
Z = Z0 = (t(Y)%*%X) #/ (n-1)
ZtZ = ZtZ0 = t(Z)%*%Z
## Total var to be explained
tvar = sum( diag ( ZtZ ) )
## Variance by PCs
r = svd(Z); d2 = r$d^2; sd2 = sum(d2); ld2 = length(d2)
var = numeric(ld2)
for (j in 1:ld2)
	var[j] = sum(d2[1:j]) / sd2
	
V = W = idx_W = NULL
bW = matrix(0, p, npcs*nlam)
cpev = numeric(npcs)

#FIXME	
## Speedup w update
#M_w = L + diag(p); res = factor_w(M_w); L_w0 = res$L; U_w0 = res$U
	
for (k in 1:npcs) {
	cat("\n\n"); print(paste("PC",k,sep=""))
	
	## Init v from SVD of data matrix
	r = svd(Z); v0 = r$v[, 1]
	y0 = Z%*%v0; Zty0 = t(Z)%*%y0; Zty0 = Zty0 / n

	#FIXME
	## Generate lambda sequence
	max_lam = max(abs(Zty0)) / alpha1
	lambdas = exp(seq(from=log(max_lam),to=log(1e-3*max_lam),length.out=nlam))
	print("lambdas")
	print(lambdas[c(1, length(lambdas))])
	
	## Search for max_lam
#	max_lam = max(abs(Zty0))
#	Cor = Zty0 / max(1, q-1) 
#	indx = which(abs(Cor)==max(abs(Cor))) 
#	for( i in seq(Cor[indx],0,length.out=50) ) {
#		print(paste("i ",i,sep=""))
#		tmp = numeric(p); tmp[indx]=i
#		lls=abs(Zty0 - ZtZ%*%tmp - alpha2*( max_lam/alpha1 )*L%*%tmp)
#		if( lls[indx] > max(lls[-indx]) ) {
#			max_lam = lls[indx] 
#			##print("t(lls)"); print(t(lls))
#			##print(paste("max_lam ", max_lam, sep=""))
#			break
#		}
#	}
#	max_lam = max_lam / alpha1
#	lambdas = exp(seq(from=log(max_lam),to=log(0.01*max_lam),length.out=nlam))
#	#lambdas = exp(seq(from=log(0.9*max_lam),to=log(1e-3*max_lam),length.out=nlam))
#	print("lambdas2")
#	print(lambdas[c(1, length(lambdas))])

	# Speedup w update
	#Z = t(YtX)
	#MZ = backsolve(U_w0, forwardsolve(L_w0, Z))
						
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
	## Method 3: Woodbury formula fast 
	#iMV = VtiMV = NULL
	#if (!is.null(V)) {
	#	myZ = Z / sqrt(n)
	#	myW = solve( diag(q) + (1/rho_v)*myZ%*%t(myZ), myZ%*%V )
	#	iMV = (1/rho_v)*V - (1/rho_v^2)*t(myZ)%*%myW
	#	VtiMV = t(iMV)%*%V
	#}
	## Method 4: Woodbury formula fast 2
	myZ = Z / sqrt(n)
	iM= (1/rho_v)*diag(p) - (1/rho_v^2)*t(myZ)%*%solve( diag(q) + (1/rho_v)*myZ%*%t(myZ) )%*%myZ
	#v
	iMV = VtiMV = NULL
	if (!is.null(V)) {
		iMV = iM%*%V
		VtiMV = t(iMV)%*%V
	}
	#w
	iMS = iM%*%S
	#StiMS = diag(ncol(S)) +  t(iMS)%*%S
	#U_w = chol(StiMS); L_w = t(U_w)
	StiMS = t(iMS)%*%S
	
	nnz_w = idx_w = 0
	#w_pls = matrix(0, p, k)
	w = u_w = u_v = numeric(p)
	vs = matrix(0, nrow=p, ncol=nlam)
	ws = matrix(0, nrow=p, ncol=nlam)
	#ys = matrix(0, nrow=q, ncol=nlam)
	
	bics = rep(1e6, nlam)
	betas = list(); prs = list(); 
	
	for (j in 1:nlam) { 

		lam = lambdas[j]; lam1 = alpha1*lam; lam2 = alpha2*lam  
		#L_w = sqrt(lam2) * L_w0; U_w = sqrt(lam2) * U_w0
		#M_w = lam2*L + rho_w*diag(p); U_w = chol(M_w); L_w = t(U_w)
		iMS = sqrt(lam2)*iMS
		D = diag(ncol(S)) + lam2*StiMS
		U_w = chol(D); L_w = t(U_w)

		# Warm start w
		v = v0; y = y0; Zty = Zty0
	
		indw = 1; iter = 1; 
		while (iter<maxit) {
			oldw = w
			#oldnnz_w = nnz_w
			#oldidx_w = idx_w
			#oldw_pls = w_pls

			## Update w
			res = update_w(n, Z, y, Zty, NULL, L, S, iS, nS, iM, iMS, L_w, U_w, oldw, NULL, rho_w, lam1, lam2) 
			x_w = res$x; w = w2 = res$z; u_w = res$u; #obj=res$obj

			if (norm(w) == 0) 
				break
			
			## Test stop condition 
			## Method 1
			indw = norm(w - oldw) / norm(oldw)
			## Method 2
			#idx_w = which(w!=0); nnz_w = length(idx_w)
        		#idx = union( idx_W, idx_w )
        		#sX = X[,idx,drop=FALSE]
        		#plsfit = pls::plsr( Y~sX, ncomp=min(k,length(idx)), method="simpls", scale=F)
			#w_pls = matrix(0, p, plsfit$ncomp); w_pls[idx, ] = plsfit$projection
			#indw = norm(w_pls - oldw_pls) / norm(oldw_pls)
			## Method 3
			#idx_w = which(w!=0); nnz_w = length(idx_w)
			#indw = (nnz_w == oldnnz_w)
			#indw =  setequal(idx_w, oldidx_w)
			
			if (indw<eps) {		
			#if (indw==1) {
				idx_w = which(w!=0); nnz_w = length(idx_w)
				
				## PLS re-estimation of w 
        			idx = union( idx_W, idx_w )
        			sX = X[,idx,drop=FALSE]
        			plsfit = pls::plsr( Y~sX, ncomp=min(k,length(idx)), method="simpls", scale=F)
        			betas[[j]] = matrix( 0, p, q ); betas[[j]][idx, ] = matrix( coef(plsfit), length(idx), q )
        			prs[[j]] = matrix( 0, p, plsfit$ncomp ); prs[[j]][idx, ] = plsfit$projection
				#prs[[j]] = w_pls

				## err term
				#beta = t(betas[[j]])
				#beta = scale(beta, F, 1/normx)
				#beta = scale(beta, -meanx, F)
				err = ( norm( Y - X%*%betas[[j]] )^2 ) / (n*q)
				
				## df term
				## Method 1
				dfterm = 0.5*(log(n*q) / (n*q)) * nnz_w
				## Method 2
        			#sX = X[,idx_w,drop=FALSE]
				#lam = lambdas[j]; lam2 = alpha2*lam
				#DF = sX %*% solve( t(sX)%*%sX + lam2*L[idx_w, idx_w] + lam2*diag(nnz_w), t(sX))
				#df =  sum(diag(DF)); dfterm = (log(n*q) / (n*q))*df
				
				## BIC
				#plty = 2*0.5*log(choose(p, nnz_w))
				bics[j] = log(err) + dfterm # + plty
				
				#DEBUG
				print(paste("[bic]j=",j,", iter=",iter,", indw=",indw,", nnz_w=",nnz_w,", err=",err,", log(err)=",log(err),", dfterm=",dfterm, ", bic=",bics[j],sep=""))
				break
			}

			## Update v
			y = Z%*%w; Zty = t(Z)%*%y; Zty = Zty / n
			res = update_v(Z, y, Zty, NULL, V, iM, iMV, VtiMV, v, NULL, rho_v)
			x_v = res$x; v = res$z; u_v = res$u
				
			## Update y
			y = Z%*%v; Zty = t(Z)%*%y; Zty = Zty / n
				
			iter = iter + 1;
		} #end while
	
		#if (DEBUG) 
			print(paste("[opt]j=",j,", iter=",iter,", indw=",indw,", nnz(w)=",length(which(w!=0)),sep=""))

		ws[ ,j] = w; vs[ ,j] = v; #ys[ ,j] = y

	} #end for 		
	
	#FIXME (shouldn't have to test this kind of condition)
	if (!length(prs)) {
		print(paste("No convergence at PC ",k,sep=""))
		break
	}
	
	idx_bic = which.min(bics)	
	w = ws[ ,idx_bic]; v = vs[ ,idx_bic]
	pr = prs[[idx_bic]]; betapls = betas[[idx_bic]]
	#DEBUG 
	print(paste("idx_bic=",idx_bic,", best nnz=",sum(w!=0),", best bic=",bics[idx_bic], sep=""))
	print("nnz(w)"); print(which(w!=0))
		
	## Update V
	V = cbind(V, v)
	## Update W  
	w = w / norm(w) 
	W = cbind(W, w)
	idx_W = union(idx_W, which(w!=0))
	## Update bW
	ws = ws %*% diag(sqrt(colSums(ws^2)))
	bW[ ,((k-1)*nlam + 1):(k*nlam)] = ws
	
	## Compute cpev
	prW = pr %*% solve( t(pr)%*%pr ) %*% t(pr) 
	cpev[k] = sum( diag( crossprod( Z0%*%prW ) ) ) / tvar
		
	#FIXME (V is orthogonal already)
	## Update data matrix
	prV = V %*% solve( t(V)%*%V ) %*% t(V) 
	Z = Z0 - Z0 %*% prV
	#ZtZ = t(Z)%*%Z
	
} #end for k
	
# Compute betapls
#TT = XO%*%W
#P = ginv(crossprod(TT)) %*% t(TT) %*% XO 
#Q = ginv(crossprod(TT)) %*% t(TT) %*% YO
#betapls = W%*%Q
	
# DEBUG
cat("\n"); print("VtV")
print(round(crossprod(V), 4))
	
## Print reg paths for pcs
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
	axis(4, at = coefs[nlam, idx], labels = colnames(X)[idx], cex = 0.1, adj=0, las=3, mar=c(8,8,1,4))
	abline(h = 0, lty = 3)
	dev.off()

	# write non-zero selected coeffs to file
	idx = which(W[,k]!=0)
	writeLines(colnames(X)[idx], paste(outdir, "pc",k,"_coeffs.txt",sep=""))
}

object <- list(X=X,Y=Y,V=V,W=W,bW=bW,betapls=betapls,var=var,cpev=cpev,meanx=meanx,normx=normx,meany=meany,normy=normy)
class(object) <- "respls"
object

}

predict.respls <- function(obj, newX) {
	
	print("[predict.respls]dim(newX)")
	print(dim(newX))

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

