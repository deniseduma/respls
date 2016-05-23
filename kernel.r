#!/usr/bin/Rscript

library(tseries)

rm(list=ls())
setwd(".")
print(getwd())

laplacian = FALSE

if (!laplacian) {
	L = NULL
	S = NULL
} else {	
	L = read.table(paste("../data/cc_regen_laplace2.csv", sep = ""), sep="\t")
	S = read.table(paste("../data/cc_regen_laplace_half2.csv", sep = ""), sep="\t")
}

if (!laplacian) 
	geneSubset = readLines("../data/inter_by_both2.txt")
if (laplacian) 
	geneSubset = readLines("../data/inter_by_both_laplace2.txt")

## read cc
regen_cc = read.table(paste("../data/cc_regen.csv", sep = ""), sep=",", header=T, row.names=1)
apoints_cc = readLines(paste("../data/cc_timepoints.txt", sep=""))
cpoints_cc = readLines(paste("../data/cc_timepoints_c.txt", sep=""))
tpoints_cc = readLines(paste("../data/cc_timepoints_t.txt", sep=""))
colnames(regen_cc) = apoints_cc

## read uu
regen_uu = read.table(paste("../data/uu_regen.csv", sep = ""), sep=",", header=T, row.names=1)
apoints_uu = readLines(paste("../data/uu_timepoints.txt", sep=""))
cpoints_uu = readLines(paste("../data/uu_timepoints_c.txt", sep=""))
tpoints_uu = readLines(paste("../data/uu_timepoints_t.txt", sep=""))

regen_cc = regen_cc[which(rownames(regen_cc) %in% geneSubset), ]
regen_uu = regen_uu[which(rownames(regen_uu) %in% geneSubset), ]

regen_big = cbind(regen_cc, regen_uu)
regen_big = t(regen_big)
regen_big = log10(1 + regen_big)

regen_t =  regen_big[which(rownames(regen_big) %in% c(tpoints_cc)), ]
regen_c =  regen_big[which(rownames(regen_big) %in% c(cpoints_cc)), ]

bs1_rows = grep("bs2", rownames(regen_t), invert=TRUE)
bs2_rows = grep("bs2", rownames(regen_t), invert=FALSE)

regen_t = (regen_t[bs1_rows, ] + regen_t[bs2_rows, ]) / 2
regen_c = (regen_c[bs1_rows, ] + regen_c[bs2_rows, ]) / 2

regen_t = scale(regen_t, center=TRUE, scale = FALSE)
regen_c = scale(regen_c, center=TRUE, scale = FALSE)

var_t = apply(regen_t, 2, function(l){var(l)})
var_c = apply(regen_c, 2, function(l){var(l)})
res_t = sort(var_t, decreasing=TRUE, index.return=TRUE)
res_c = sort(var_c, decreasing=TRUE, index.return=TRUE)
lim = min(length(geneSubset), 2000)
regen_t = regen_t[, res_t$ix[1:lim]]
regen_c = regen_c[, res_t$ix[1:lim]]

timepts = rownames(regen_t)
genes = colnames(regen_t)

dim = dim(regen_t)
n = dim[1]
p = dim[2]
print("dim(regen_t)")
print(dim(regen_t))

#for (i in 1:5) {
#	g <- regen_t[, i]
#	r.arma <- arma(g, order = c(1, 0), include.intercept = T)
#	print(summary(r.arma))
	#print(coef(r.arma))
	#plot(r.arma)
#}

#Estimate AR(1) model for each gene
sigma2 <- 0
a <- numeric(p)
for (i in 1:p) {
	num <- 0
	denom <- 0
	for (t in 2:n) {
		num <- num + regen_t[t-1, i] * regen_t[t, i]
		denom <- denom + regen_t[t-1, i] * regen_t[t-1, i]
	}
	a[i] <- num/denom
	
	s <- 0
	for (t in 2:n) {
		e <- regen_t[t, i] - a[i] * regen_t[t-1, i]
		s <- s + e * e
	}
	s <- s / (n-1)
	sigma2 <- sigma2 + s
}
sigma2 <- sigma2 / p

#Compute kernel matrix
lambda <- 0.5
expl <- exp(-lambda)
K <- matrix(0, p, p)
for (i in 1:p) {
	for (j in 1:p) {
		num <- expl * a[i] * a[j]
		denom <- 1 - num
		m <- num / denom
		m2 <- 1 / denom
		K[i,j] <- regen_t[1, i] * m * regen_t[1, j] + 1/(exp(lambda) - 1) * sigma2 * m2	
	}
}

#Write kernel matrix to file
write.table(K, file = "../data/inter_by_both2_kernel.txt", quote = F, sep = " ", row.names = F, col.names = F)

