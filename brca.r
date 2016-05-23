#!/usr/bin/Rscript


rm(list=ls())
setwd(".")
print(getwd())

library(parallel)

source("nrm.r")
source("opt.r")
source("ust.R")
source("spls.R")
source("correctp.R")
source("spls.dv.R")
source("predict.spls.R")

indir = "brca/"
outdir = "brca/"

method = "respls" 


##Read data
X = as.matrix(read.table(paste(indir, "brca.expr.dat3.txt",sep=""), header=T, row.names=1, check.names=F, stringsAsFactors=F, sep=" "))
Y = as.matrix(read.table(paste(indir, "rppa.dat2.txt", sep=""), header=T, row.names=1, check.names=F, stringsAsFactors=F, sep=" "))
X = log(1+X)
Y = Y[ ,1:10]
print("str(X)")
print(str(X))
print("str(Y)")
print(str(Y))

n = nrow(X)
p = ncol(X)
q = ncol(Y)
one = rep(1, n)

## Check normality of X and Y
#dev.new()
#X = scale(X,T,F)
#hist(X, main="X, mean 0")
#dev.new()
#hist(Y, main="Y")

npcs=10;n=nrow(X);p=ncol(X);q=ncol(Y)

## Adj matrix
A =  as.matrix(read.table(paste(indir,"adj.txt",sep=""), header = T, row.names=1, check.names=T, stringsAsFactors=F, sep=" "))
print("dim(A)")
print(dim(A))

## Proper Laplacian
D <- diag(rowSums(A))
L <- D - A
#d <- diag(D)
#d[d>0] <- 1 / sqrt(d[d>0])
#L <- diag(d) %*% L %*% diag(d)
## Placeholder Laplacian
#L=diag(p)

nlam = 100
##Run spls/respls 
start <- proc.time()
if (method == "spls") {
	res = spls(X, Y, K=npcs, eta=0.5, kappa=0.5, select="pls2", fit="simpls", scale.x=TRUE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=FALSE)
} else if (method == "respls") {
	res = opt(X, Y, L, NULL, NULL, 0, npcs, nlam,outdir)
}
#Measure execution time
elapsed <- proc.time() - start
print(paste("execution time ", round(elapsed[3], 4), sep=""))
#print("res")
#print(res)

if (method == "spls") {
	#betapls <- res$betahat
	#W = sapply(res$ws, unlist)
	W = as.matrix(res$ws[[1]])
	if (length(res$ws)>1)
		for (j in 2:length(res$ws)) 
			W = cbind(W, res$ws[[j]])
	# Normalize direction vectors 
	W = W %*% diag( 1/sqrt(colSums(W^2)), nrow=ncol(W), ncol=ncol(W) ) 
} else if (method == "respls") {
	W = res$W
}

## cpev
if (method=="respls") {
	print("var")
	print(res$var)
	print("cpev")
	print(res$cpev)
}

## Print gene list to file
#for (j in 1:ncol(W)) {
#	writeLines(colnames(X)[which(W[ ,j]!=0)], paste(outdir,method, "_pc",j,"_coeffs.txt",sep="")) 
#}

