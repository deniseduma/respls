#!/usr/bin/Rscript


rm(list=ls())
setwd(".")
print(getwd())

library(multicore)

source("nrm.r")
source("runParallel1.r")

method = "respls" #"spls"

outdir = "sim1/"

npcs = 1; n=40; p=1000; q=1; yprop = 0.5
mu1 = c(rep(-2, n/2), rep(2, n/2))
mu2 = c(rep(-1, n/4), rep(1, n/4), rep(-1, n/4), rep(1, n/4))
p1 = 200; p2 = 50

# Ground truth
v1 = c(rep(1,p1), rep(0, p-p1))
v2 = c(rep(0, p1), rep(1, p2), rep(0, p-p1-p2))
v3 = c(rep(1,p1+p2), rep(0, p-p1-p2))
v1 = v1/norm(v1) 
v2 = v2/norm(v2)
v3 = v3/norm(v3)
## yprop = 0
#V0 = cbind(v2)
## yprop = 1
#V0 = cbind(v1)
## yprop = 0.5
V0 = cbind(v3)

# Build Laplacian on genes 201-250
A <- matrix(0, p, p)
## yprop = 0
#A[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- 1 
## yprop = 1
#A[1:p1, 1:p1] <- 1 
## yprop = 0.5
A[1:p1, 1:p1] <- 1 
A[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- 1 
diag(A) <- 0
D <- diag(rowSums(A))
L <- D - A
#d <- diag(D)
#d[d>0] <- 1/sqrt(d[d>0])
#L <- diag(d) %*% L %*% diag(d)
## Placeholder Laplacian
L=diag(p)

nlam = 100
num_iter = 50
res = mclapply(1:num_iter, function(i) runParallel1(i,n,p,p1,p2,q,yprop,mu1,mu2,npcs,nlam,V0,L,outdir,method), mc.preschedule=FALSE, mc.cores = 15)
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
prec = prec / num_iter;
recall = recall / num_iter;

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
	elapsed = elapsed + crt$elapsed
}
elapsed = elapsed / num_iter
cat("\n")
print(paste("execution time ", round(elapsed[3], 4), sep=""))


# csplit
csplit = 0
for (i in 1:num_iter) {
	crt = res[[i]]
	csplit = csplit + crt$csplit
}
cat("\n")
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
print("sds_recall")
print(sds_recall)
cat("\n")
print("avgs_prec")
print(avgs_prec)
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
#	plot(prec[j, ], recall[j, ], cex=1, pch=16, type="l", xlim=c(0, 0.001), ylim=c(0,1), xlab ="FPR", ylab="TPR", main=paste("ROC curve PC", j,sep=""), cex.lab=1, cex.axis=1)
#	points(prec[j, ], recall[j, ])
#	dev.off() 
#}
#}

