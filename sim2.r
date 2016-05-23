#!/usr/bin/Rscript


rm(list=ls())
setwd(".")
print(getwd())

library(multicore)

source("nrm.r")
source("runParallel2.r")

method = "respls" #"spls"

outdir = "sim2/"

npcs = 2; n=20; p=100; q=2
p1 = 0.4*p; p2 = 0.8*p
c1 = 290; c2 = 300; c3 = 0.3*0.3*c1 + 0.925*0.925*c2

#Compute covariance matrix
C = matrix(0, p, p)
C[1:p1, 1:p1] = c1
C[(p1+1):p2, (p1+1):p2] = c2
C[(p2+1):p, (p2+1):p] = c3
C[1:p1, (p2+1):p] = -0.3*c1
C[(p2+1):p, 1:p1] = -0.3*c1
C[(p1+1):p2, (p2+1):p] = 0.925*c2
C[(p2+1):p, (p1+1):p2] = 0.925*c2
C = C + diag(p)

#Useful to compute tpr and fpr
v2 = c(rep(1,p1), rep(0,p-p1))
#v1 = c(rep(0,p1), rep(1,p2-p1), rep(0,p-p2))
v1 = c(rep(0,p1), rep(1,p-p1))
v3 = c(rep(0,p2), rep(1,p-p2))
v1 = v1/norm(v1) 
v2 = v2/norm(v2)
v3 = v3/norm(v3)
V0 = cbind(v1, v2)

## Adjacency matrix
#A <- cov2cor(C)
A <- abs(solve(C))
diag(A) <- 0
## Proper Laplacian
D <- diag(rowSums(A))
L <- D - A
#DD <- diag(diag(D)^(-1/2))
#L <- DD %*% L %*% DD
## Placeholder Laplacian
L=diag(p)

nlam = 100
num_iter = 50
res = mclapply(1:num_iter, function(i) runParallel2(i,n,p,p1,p2,q,c1,c2,c3,npcs,nlam,C,V0,L,outdir,method), mc.preschedule=FALSE, mc.cores = 15)
print("str(res)")
print(str(res))

prec = matrix(0, npcs, nlam)
recall = matrix(0, npcs, nlam)
prec2 = matrix(0, npcs, num_iter)
recall2 = matrix(0, npcs, num_iter)


if (method=="respls") {
for (i in 1:num_iter) {
	crt = res[[i]]
	for (j in 1:npcs) { 
		for (l in 1:nlam) {
			prec[j,l] = prec[j,l] + crt$prec[j, l]
			recall[j,l] = recall[j,l] + crt$recall[j, l]
		}
	}
}	
}
prec = prec / num_iter
recall = recall / num_iter


for (i in 1:num_iter) {
	crt = res[[i]]
	for (j in 1:npcs) { 
		prec2[j,i] = crt$prec2[j]
		recall2[j,i] = crt$recall2[j]
	}
}

# tvar
if (method=="respls") {
tvar = matrix(0, num_iter, q)
for (i in 1:num_iter) {
	crt = res[[i]]
	tvar[i, ] = crt$var
}
avgs_tvar = numeric(q)
sds_tvar = numeric(q)
for (j in 1:q) { 
	avgs_tvar[j] = mean(tvar[ ,j])
	sds_tvar[j] = sd(tvar[ ,j])
}
cat("\n")
print("avgs_tvar")
print(avgs_tvar)
print("sds_tvar")
print(sds_tvar)

# cpev
nlen = length(res[[1]]$cpev)
cpev = matrix(0, num_iter, nlen)
for (i in 1:num_iter) {
	crt = res[[i]]
	cpev[i, ] = crt$cpev
}	
avgs_cpev = numeric(nlen)
sds_cpev = numeric(nlen)
for (j in 1:nlen) { 
	avgs_cpev[j] = mean(cpev[ ,j])
	sds_cpev[j] = sd(cpev[ ,j])
}
print("avgs_cpev")
print(avgs_cpev)
print("sds_cpev")
print(sds_cpev)
}

# MSPE
mspe = matrix(0, num_iter, q)
for (i in 1:num_iter) {
	crt = res[[i]]
	mspe[i, ] = crt$mspe
}
avgs_mspe = numeric(q)
sds_mspe = numeric(q)
for (j in 1:q) { 
	avgs_mspe[j] = mean(mspe[ ,j])
	sds_mspe[j] = sd(mspe[ ,j])
}
cat("\n")
print("avgs_mspe")
print(avgs_mspe)
print("sds_mspe")
print(sds_mspe)

# elapsed
elapsed = 0
for (i in 1:num_iter) {
	crt = res[[i]]
	elapsed <- elapsed + crt$elapsed
}
elapsed = elapsed/num_iter
cat("\n")
print(paste("execution time ", round(elapsed[3], 4), sep=""))


# csplit
cat("\n")
csplit = 0
for (i in 1:num_iter) {
	crt = res[[i]]
	csplit = csplit + crt$csplit
}
csplit = csplit / num_iter
print(paste("csplit ",csplit,sep=""))

# t_prec
t_prec = numeric(num_iter)
for (i in 1:num_iter) { 
	crt = res[[i]]
	t_prec[i] = crt$t_prec
}	
avgs_tprec = mean(t_prec)
sds_tprec = sd(t_prec)
cat("\n")
print("avgs_tprec")
print(avgs_tprec)
print("sds_tprec")
print(sds_tprec)


#if (method=="respls") {
#cat("\n")
#print("fpr and tpr")
#for (j in 1:npcs) { 
	#print(paste("roc_fpr[", j, "]", sep=""))
	#print(roc_fpr[j, ])
	#print(paste("roc_tpr[", j, "]", sep=""))
	#print(roc_tpr[j, ])
#}
#}

# recall and prec
avgs_recall = numeric(npcs)
avgs_prec = numeric(npcs)
sds_recall = numeric(npcs)
sds_prec = numeric(npcs)
for (j in 1:npcs) { 
	avgs_recall[j] = mean(recall2[j,])
	avgs_prec[j] = mean(prec2[j,])
	sds_recall[j] = sd(recall2[j,])
	sds_prec[j] = sd(prec2[j,])
}

cat("\n")
print("avgs_recall")
print(avgs_recall)
print("avgs_prec")
print(avgs_prec)
print("sds_recall")
print(sds_recall)
print("sds_prec")
print(sds_prec)

#ma = matrix(0,npcs,num_iter)
#for (i in 1:num_iter) {
#	crt = res[[i]]
#	for (j in 1:npcs) 
#		ma[j,i] = crt$ma_angle[j]
#}		
#cat("\n")
#for (j in 1:npcs) {
#	print(paste("Median angle for pc ", j, sep=""))
#	print(median(ma[j, ]))
#}

#Plot ROC curves for each PC
#if (method=="respls") {
#for (j in 1:npcs) { 
#	pdf(paste("pc", j,"_roc.pdf",sep=""))
#	plot(roc_fpr[j, ], roc_tpr[j, ], cex=1, pch=16, type="l", xlim=c(0, 0.001), ylim=c(0,1), xlab ="FPR", ylab="TPR", main=paste("ROC curve PC", j,sep=""), cex.lab=1, cex.axis=1)
#	points(roc_fpr[j, ], roc_tpr[j, ])
#dev.off() 
#}
#}

