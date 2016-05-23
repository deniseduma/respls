#!/usr/bin/Rscript
rm(list=ls())
setwd(".")
print(getwd())

library(multicore)

source("nrm.r")
source("runParallelSPLS.r")
source("runParallelRESPLS.r")

method = "respls" #"spls"

outdir = "rex1/"

split = 40

# Read X 
X = read.csv("rex1/xmatrix.csv", header = TRUE, row.names = 1, sep = ",", quote = "")
X = as.matrix(t(X))
print("dim(X)")
print(dim(X))

# Read Y
Y = read.csv("rex1/ymatrix.csv", header = TRUE, row.names = 1, sep = ",", quote = "\"")
rnames = rownames(Y)
cnames = colnames(Y)
Y = as.matrix(Y[, "cag"])
rownames(Y) = rnames
colnames(Y) = c("cag")
print("dim(Y)")
print(dim(Y))
print(head(Y, 5))

# Build Laplacian
A <- readRDS("rex1/x_lapl_adj.rds")
D <- diag(rowSums(A))
print("dim(A)")
print(dim(A))
L <- D - A
#d <- diag(D)
#d[d>0] <- 1/sqrt(d[d>0])
#L <- diag(d) %*% L %*% diag(d)
# Placeholder Laplacian
L <- diag(ncol(A))

X <- X[, which(colnames(X) %in% colnames(A))]
print("dim(X)")
print(dim(X))

# Retain top variable genes
X.var <- apply(X, 2, var)
idx_x = ( X.var >= quantile(X.var, 0.5, na.rm=T)  & !is.na(X.var))
X = X[ ,idx_x]
L = L[idx_x, idx_x]
print("dim(X)")
print(dim(X))
print("dim(L)")
print(dim(L))

npcs = 5; n = nrow(X); p = ncol(X); q = ncol(Y)  

nlam = 100
num_iter = 1
if (method=="spls") {
	res <- mclapply(1:num_iter, function(i) runParallelSPLS(i,X,Y,npcs,split), mc.preschedule=FALSE, mc.cores = 15)
} else if (method=="respls") {
	res <- mclapply(1:num_iter, function(i) runParallelRESPLS(i,X,Y,npcs,nlam,L,split), mc.preschedule=FALSE, mc.cores = 15)
}
print("str(res)")
print(str(res))

#res=res[-(which(sapply(res,is.null),arr.ind=TRUE))]
res = res[lapply(res, length) != 0]
if (length(res)==0)
	stop("No results!")
num_iter = length(res)
print("length(res)")
print(length(res))

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
cpev = matrix(0, num_iter, npcs)
for (i in 1:num_iter) {
	crt = res[[i]]
	cpev[i, ] = crt$cpev
}	
avgs_cpev = numeric(npcs)
sds_cpev = numeric(npcs)
for (j in 1:npcs) { 
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


#Plot ROC curves for each PC
#if (method=="respls") {
#for (j in 1:npcs) { 
#	pdf(paste("pc", j,"_roc.pdf",sep=""))
#	plot(prec[j, ], recall[j, ], cex=1, pch=16, type="l", xlim=c(0, 0.001), ylim=c(0,1), xlab ="FPR", ylab="TPR", main=paste("ROC curve PC", j,sep=""), cex.lab=1, cex.axis=1)
#	points(prec[j, ], recall[j, ])
#	dev.off() 
#}
#}



