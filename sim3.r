#!/usr/bin/Rscript


rm(list=ls())
setwd(".")
print(getwd())

library(Matrix)
library(parallel)

source("nrm.r")
source("runParallel3.r")

indir = "sim3/"
outdir = "sim3/"

method = "respls" 

sim4 = 1
sim5 = 0

npcs=5;n=50;p=500;pp1=10;pp2=100;q=5

## The eigenvalues of the covariance matrix
cc = c(600, 500, 400, 300, 200, rep(1, (p-npcs))) 
if (sim4) {
	cc =  read.table(paste(indir,"c_sizes.txt",sep=""), header = FALSE, quote = NULL)
	cc = as.numeric(cc[, 1])
	npcs=q=length(cc) 
	#npcs=q=2 
	cc = c(cc[1:npcs], rep(1, p-npcs))
	print("cc"); print(cc)
	pcs = seq(1:10)
}

## The eigenvectors of the covariance matrix
v1 = c(rep(1, pp1), rep(0, (p-pp1)))
v2 = c(rep(0, pp2), rep(1, pp1), rep(0, (p-pp2-pp1)))
v3 = c(rep(0, (2*pp2)), rep(1, pp1), rep(0, (p-2*pp2-pp1)))
v4 = c(rep(0, (3*pp2)), rep(1, pp1), rep(0, (p-3*pp2-pp1)))
v5 = c(rep(0, (4*pp2)), rep(1, pp1), rep(0, (p-4*pp2-pp1)))
#v1 = v1/sqrt(sum(v1^2)) 
#v2 = v2/sqrt(sum(v2^2))
#v3 = v3/sqrt(sum(v3^2))
#v4 = v4/sqrt(sum(v4^2))
#v5 = v5/sqrt(sum(v5^2))
V0 = cbind(v1,v2,v3,v4,v5) 
if (sim4) {
	V0 =  as.matrix(read.table(paste(indir,"V.txt",sep=""), header = FALSE, quote = NULL))
	V0 = V0[ ,1:npcs, drop=F]
	cat("\n"); print("dim(V0)"); print(dim(V0))
	for ( j in 1:npcs) {
		idx = which(V0[ ,j]>0)
		print(paste("[PC",j,"], nnz=",length(idx),sep=""))
		print(idx)
	}
}

## Remaining eigenvectors
V = V0
for (i in (npcs+1):p) {
	#v = rnorm(p)
	v = runif(p, 0, 1)
	V = cbind(V, v)
}
VtV = t(V) %*% V
#print("VtV")
#print(diag(VtV))
#print(VtV[1:100,1:100])
V = qr.Q(qr(V))

# Covariance matrix
C = matrix(0,p,p)
for (i in 1:p) {
	C = C + cc[i] * V[, i] %*% t(V[, i])
}

## Adj matrix
#A <- matrix(0, p, p)
#A[1:pp1, 1:pp1] <- 1 
#A[(pp2+1):(pp2+pp1), (pp2+1):(pp2+pp1)] <- 1 
#A[(2*pp2+1):(2*pp2+pp1), (2*pp2+1):(2*pp2+pp1)] <- 1 
#A[(3*pp2+1):(3*pp2+pp1), (3*pp2+1):(3*pp2+pp1)] <- 1 
#A[(4*pp2+1):(4*pp2+pp1), (4*pp2+1):(4*pp2+pp1)] <- 1 
#diag(A) <- 0
if (sim4) {
	A =  as.matrix(read.table(paste(indir,"A.txt",sep=""), header = FALSE, quote = NULL))
	print("dim(A)")
	print(dim(A))
}

## Proper Laplacian
D <- diag(rowSums(A))
L <- D - A
#d <- diag(D)
#d[d>0] <- 1/sqrt(d[d>0])
#L <- diag(d) %*% L %*% diag(d)
#if (sim4) {
#	L =  as.matrix(read.table(paste(indir,"L.txt",sep=""), header = FALSE, quote = NULL))
#	print("dim(L)")
#	print(dim(L))
#}
## Placeholder Laplacian
L=diag(p)

r = eigen(L)
S = r$vectors%*%diag ( sqrt(r$values) )
iS = diag ( 1/sqrt(r$values))  %*% t(r$vectors)
nS = norm(S)^2
#print("class(S)")
#print(class(S))
#print("dim(S)")
#print(dim(S))
#print("rank(L)")
#print(rankMatrix(L))
#print("rank(S)")
#print(rankMatrix(S))

nlam = 50
num_iter = 1
if (sim5) {
	pcs = c(1,  3, 5,  7)
	q = length(pcs)
}

## Run spls/respls 
#try(
	#res = mclapply(1:num_iter, function(i) runParallel3(i,n,p,q,npcs,nlam,C,V0,L,S,iS,nS,pcs,outdir,method), mc.preschedule=FALSE, mc.cores = 15)
#)
	res = lapply(1:num_iter, function(i) runParallel3(i,n,p,q,npcs,nlam,C,V0,L,S,iS,nS,pcs,outdir,method))
#print("res")
#print(res)

X = res[[1]]$X
Y = res[[1]]$Y
U = res[[1]]$U
#if (method=="respls") {
#	P = res[[1]]$P
#	Q = res[[1]]$Q
#}

#res=res[-(which(sapply(res,is.null),arr.ind=TRUE))]
#res = res[lapply(res, length) != 0]
#if (length(res)==0)
#	stop("No results!")
#num_iter = length(res)
#print("length(res)")
#print(length(res))

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

## cpev
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
#csplit = 0
#for (i in 1:num_iter) {
#	crt = res[[i]]
#	csplit = csplit + crt$csplit
#}
#cat("\n")
#csplit = csplit / num_iter
#print(paste("csplit ",csplit,sep=""))

# t_prec
#t_prec = numeric(num_iter)
#for (i in 1:num_iter) { 
#	crt = res[[i]]
#	t_prec[i] = crt$t_prec
#}	
#avgs_tprec = mean(t_prec)
#sds_tprec = sd(t_prec)
#cat("\n")
#print("avgs_tprec")
#print(avgs_tprec)
#print("sds_tprec")
#print(sds_tprec)

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

