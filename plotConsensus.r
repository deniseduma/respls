#!/usr/bin/Rscript

rm(list=ls())
setwd(".")

library(gplots)

myncomp = 8
tdata = read.table(paste("../data/cc/consensus_k", myncomp, "_t_3.txt", sep=""), header = FALSE, sep = " ", quote = NULL)
cdata = read.table(paste("../data/cc/consensus_k", myncomp, "_c_3.txt", sep=""), header = FALSE, sep = " ", quote = NULL)
print(dim(tdata))
print(dim(cdata))

tdata[tdata>1] = 1
cdata[cdata>1] = 1

##Treated
sum1 = apply(tdata, 1, sum)
idx1 = which(sum1 == 0)
tdata = tdata[-idx1, -idx1]
cat("\ndim(tdata) after dim red\n")
cat(dim(tdata))

##Control
sum1 = apply(cdata, 1, sum)
idx1 = which(sum1 == 0)
cdata = cdata[-idx1, -idx1]
cat("\ndim(cdata) after dim red\n")
cat(dim(cdata))

options("expressions"=5e4)

myDist = function(x) {as.dist(1-x)}
myHclustDist = function(x) { hclust(d=x, method="average") }

pdf(paste("../data/cc/consensus_k", myncomp,"_t_3_max1.pdf",sep=""))
heatmap.2(as.matrix(tdata), dendrogram="column", symm=TRUE, dist = myDist, hclust = myHclustDist, col=greenred(75), key=TRUE, density.info="none", trace="none", main=paste("Treated consensus matrix, k=", myncomp,sep=""))
dev.off()

pdf(paste("../data/cc/consensus_k", myncomp,"_c_3_max1.pdf",sep=""))
heatmap.2(as.matrix(cdata), dendrogram="column", symm=TRUE, dist = myDist, hclust = myHclustDist, col=greenred(75), key=TRUE, density.info="none", trace="none", main=paste("Control consensus matrix, k=", myncomp,sep=""))
dev.off()
