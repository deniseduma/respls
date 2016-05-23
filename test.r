library(MASS)
source("update_w.r")
library(glmnet)

n=2; p=100; q=1
X <- mvrnorm(n, rep(0, p), diag(p))
Y <- rnorm(n)
X = scale(X,T,T)
Y = scale(Y,T,T)

Lapl = diag(p)

alpha = 0.8
ac = (t(X)%*%Y) / (n-1)
indx = which(abs(ac)==max(abs(ac)))
print("ac")
print(t(ac))
print(paste("max(abs(ac)) ",ac[indx],", indx ",indx,sep=""))
#ac = cor(X, Y)
#indx = which(abs(ac)==max(abs(ac))) #asign = sign(ac[indx])
#print(paste("indx ",indx,", ac2 ",ac[indx], sep=""))
for( i in seq(ac[indx],0,length.out=50) ) {
	print(paste("i ",i,sep=""))
	tmp = numeric(p); tmp[indx]=i
	lls=abs(t(X)%*%(Y-X%*%tmp)-(1-alpha)*( max(abs(t(X)%*%Y)) / alpha )*Lapl%*%tmp)
	if(lls[indx] > max(lls[-indx])) {
		max_l = lls[indx] 
		#print("t(lls)"); print(t(lls))
		print(paste("max_l1 ", max_l, sep=""))
		break;
	}
}
nlam=5; max_l = max_l/alpha
print(paste("max_l2 ", max_l, sep=""))
#max_l = max(abs(t(X)%*%Y))/alpha
lambdas = exp(seq(from=log(max_l),to=log(0.01*max_l),length.out=nlam))
print("lambdas"); print(t(lambdas))

start = proc.time()

Z = NULL
###U = chol( diag(p) + Lapl ); L = t(U)
#U = chol( Lapl ); L = t(U)
for (j in 1:nlam) {
	lam1 = alpha*lambdas[j]; lam2 = (1-alpha)*lambdas[j]
	#L = sqrt(lam2)*L; U = sqrt(lam2)*U
	U = chol( diag(1,p,p) + lam2*Lapl ); L = t(U)
	res = update_w(A=X,b=Y,Atb=t(X)%*%Y,AtA=t(X)%*%X,Lapl=Lapl,S=NULL,invS=NULL,nrmS=0,L,U,Zhat=NULL,D=NULL,z_init=numeric(p),u_init=numeric(p),rho=1,lam1=lam1,lam2=lam2)
	z = res$z; Z = cbind(Z, z)
}

elapsed = proc.time() - start
print( paste("elapsed ",elapsed[3], sep=""))

print("admm")
colnames(Z)=NULL; rownames(Z)=NULL
print(colSums(Z!=0))
apply(Z,2,function(x){print(which(x!=0))})
	
res2 = glmnet(X,Y,family="mgaussian",alpha=alpha,intercept=F,standardize=F,lambda=lambdas/n)
print("glmnet")
m=as.matrix(res2$beta);colnames(m)=NULL;rownames(m)=NULL
print(colSums(m!=0))
apply(m,2,function(x){print(which(x!=0))})

