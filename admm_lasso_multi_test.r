library("glmnet")
library("elasticnet")
source("admm_lasso_multi2.r")

m <- 200
K <- 20
ni <- 1000
n <- ni*K
p <- 10/K

rho <- 1
alpha <- 1.8

#generate block sparse solution vector
x <- matrix(0, ni , K)
for (i in 1:K) 
	if (runif(1) < p)
		x[, i] <- rnorm(ni)
xtrue <- x <- Matrix(as.vector(x))

#generat random data matrix
A <- Matrix(rnorm(m*n), nrow=m, ncol=n)
#normalize columns of A
A <- A %*% Diagonal(n, 1/sqrt(colSums(A^2)))
#generate measurement b with noise
b <- A %*% x + Matrix(rnorm(m, 0, sqrt(1)))

#lambda max
nrms = numeric(K)
for (i in 1:K) {
	Ai <- A[ ,((i-1)*ni+1):(i*ni)]
	Aitb <- t(Ai) %*% b
	nrms[i] <-  sqrt(sum(Aitb^2))
}
lam_max <- max(nrms)
#lam_max <- max(abs(t(A) %*% b))

#regularization parameter
lam <- 0.01*lam_max #best 0.01 for both admm_lasso_multi and admm_lasso_multi2

#Precompute matrix decomps
decomps <- vector("list", K)
for (i in 1:K) {
	#Ai <- A
	Ai <- A[ ,((i-1)*ni+1):(i*ni)]
	res<-factor(Ai, rho) 
	decomps[[i]]<-list(L=res$L, U=res$U)
}

init_x <- Matrix(0, n , 1)
start <- proc.time()
#x <- admm_lasso_multi(A,b,init_x,decomps,lam,rho,alpha,F,ni,20)
x <- admm_lasso_multi2(A,b,init_x,lam,rho,alpha,F,ni,20)
elapsed <- proc.time() - start
print(paste("Elapsed time ", elapsed[3], sep=""))

cat("\n")
print(paste("norm(xtrue)=", sqrt(sum(xtrue^2)), ", norm(x)=", sqrt(sum(x^2)), ", length(x)=",length(x), sep=""))

print("glmnet solution")
r <- glmnet(A, b, family="gaussian", nlambda=1, lambda=lam, standardize=FALSE, intercept=FALSE)
print(paste("norm(beta)=", sqrt(sum(r$beta^2)),", length(r$beta)=",length(r$beta),sep=""))

#print("enet solution")
#r <- enet(as.matrix(A), as.matrix(b), normalize=FALSE, intercept=FALSE)
#t <- t(r$beta.pure)
##t <- (t!=0) 
##beta <- rowMeans(t)
##beta[beta<0.6] <- 0
#print(paste("norm(beta)=", sqrt(sum(t[, ncol(t)]^2)),", length(beta)=",length(t[ , ncol(t)]),sep=""))


