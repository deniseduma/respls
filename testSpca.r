library("MASS")
library("Matrix")
library("elasticnet")

n=50; p=100;
pp=p
p1=0.4*p;
p2=0.8*p;

ncomp = 2 

c1 = 290
c2 = 300
c3 = 0.3*0.3*c1 + 0.925*0.925*c2

C = matrix(0, p, p)
C[1:p1, 1:p1] = c1
C[(p1+1):p2, (p1+1):p2] = c2
C[(p2+1):pp, (p2+1):pp] = c3
C[1:p1, (p2+1):pp] = -0.3*c1
C[(p2+1):pp, 1:p1] = -0.3*c1
C[(p1+1):p2, (p2+1):pp] = 0.925*c2
C[(p2+1):pp, (p1+1):p2] = 0.925*c2
C = C + diag(p)

X <- mvrnorm(n, rep(0, p), C) 
X <- X %*% diag(1/sqrt(colSums(X^2)))
print("dim(X)")
print(dim(X))
print(paste("cond(X)=",kappa(X,exact=TRUE),sep=""))
print(paste("rank(X)=",rankMatrix(X),sep=""))

v1 = c(rep(0,p1), rep(1,pp-p1))
v2 = c(rep(1,p1), rep(0,pp-p1))
v1 = v1/sqrt(sum(v1^2)) 
v2 = v2/sqrt(sum(v2^2))

V0 = cbind(v1, v2)
#print('V0')
#print(V0)

lam1 = matrix(c(0.01,0.01,0.05,0.05,0.1,0.1,1,1,10,10,100,100),2,6)
print("lam1")
print(lam1)

num_iter = 10

#tp = numeric(ncomp)
#fp = numeric(ncomp)
#tn = numeric(ncomp)
#fn = numeric(ncomp)
tpr = numeric(ncomp)
fpr = numeric(ncomp)
bic = numeric(ncomp)

for  (i in 1:num_iter) {
	
	min_lam1 = numeric(ncomp)
	min_ang = rep(1e3, ncomp)
	V_hat = matrix(0,p,ncomp)
	
	#try multiple lambda values
	for (j in 1:ncol(lam1)) {
		
		cat("\n")
		print(paste("j=",j,", lam1=",lam1[,j][1],sep=""))
		res=spca(X,2,para=lam1[,j],type="predictor",sparse="penalty")
	
		print("res$loadings")
		print(res$loadings)

		V_hat0 = res$loadings
		print("dim(V_hat0)")
		print(dim(V_hat0))

		#find best recovery for each PC
		for (k1 in 1:ncomp) {
			#print(paste("k1=",k1,sep=""))
			for (k2 in 1:ncomp) {
				nrm = sqrt(sum(V_hat0[,k2]^2))
				if (nrm==0)
					nrm = 1
				ang = acos(abs(crossprod(V0[,k1], V_hat0[,k2]) / nrm))
				#print(paste("k2=",k2,", ang=",ang,sep=""))
				if (ang<min_ang[k1]) {
					min_ang[k1] = ang
					V_hat[,k1] = V_hat0[,k2]
					min_lam1[k1] = lam1[,j][k1]
				}
			}
		}	
	
	}#end j
	
	cat("\n")
	print(paste("iter=",i,sep=""))
	for (k in 1:ncomp) 
	print(paste("pc",k,", best recovery ",min_lam1[k],sep=""))
	#print("V_hat")
	#print(V_hat)

	#compute tpr and fpr	
	for (k in 1:ncomp) {
		tp = sum(V0[,k] & V_hat[,k])
		fn = sum(V0[,k]!=0) - tp
		fp = sum((!V0[,k]) & V_hat[,k])
		tn = sum(!V0[,k]) - fp
		crt_fpr = fp / (fp+tn)
		crt_tpr = tp / (tp+fn)
		fpr[k] = fpr[k] + crt_fpr
		tpr[k] = tpr[k] + crt_tpr
		#print(paste("k=",k,", tp=",tp,", fp=",fp,", tn=",tn,", fn=",fn,sep=""))
		print(paste("k=",k,", tpr=",crt_tpr,", fpr=",crt_fpr,sep=""))
	}
	
}#end i

cat("\n")
tpr = tpr/num_iter
fpr = fpr/num_iter
for (k in 1:ncomp)
	print(paste("k=",k,", tpr=",tpr[k],", fpr=",fpr[k],sep=""))


