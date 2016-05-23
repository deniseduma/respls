#!/usr/bin/Rscript

library(gplots)

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
colnames(regen_uu) = apoints_uu

regen_cc = as.matrix(regen_cc)
regen_uu = as.matrix(regen_uu)

geneSubset = readLines("../data/inter_by_both.txt")
regen_cc = regen_cc[geneSubset, ]
regen_uu = regen_uu[geneSubset, ]

regen_big = cbind(regen_cc, regen_uu)
regen_big = t(regen_big)
print("dim(regen_big)")
print(dim(regen_big))

genes = colnames(regen_big)

regen_t =  regen_big[which(rownames(regen_big) %in% c(tpoints_uu)), ]
regen_c =  regen_big[which(rownames(regen_big) %in% c(cpoints_uu)), ]

regen_t = scale(regen_t, center=TRUE, scale = TRUE)
regen_c = scale(regen_c, center=TRUE, scale = TRUE)

##cc
#tpoints_t_1 = which(rownames(regen_t) %in% c("cc_0t", "cc_24t", "cc_48t", "cc_72t", "cc_96t", "cc_120t", "cc_144t", "cc_168t"))
#tpoints_c_1 = which(rownames(regen_c) %in% c("cc_0c", "cc_24c", "cc_48c", "cc_72c", "cc_96c", "cc_120c", "cc_144c", "cc_168c"))
#tpoints_t_2 = which(rownames(regen_t) %in% c("cc_0t_bs2", "cc_24t_bs2", "cc_48t_bs2", "cc_72t_bs2", "cc_96t_bs2", "cc_120t_bs2", "cc_144t_bs2", "cc_168t_bs2"))
#tpoints_c_2 = which(rownames(regen_c) %in% c("cc_0c_bs2", "cc_24c_bs2", "cc_48c_bs2", "cc_72c_bs2", "cc_96c_bs2", "cc_120c_bs2", "cc_144c_bs2", "cc_168c_bs2"))
##uu
tpoints_t_1 = which(rownames(regen_t) %in% c("0t", "24t", "48t", "54t", "60t", "66t", "72t", "96t", "120t", "144t", "168t"))
tpoints_c_1 = which(rownames(regen_c) %in% c("0c", "24c", "48c", "54c", "60c", "66c", "72c", "96c", "120c", "144c", "168c"))
tpoints_t_2 = which(rownames(regen_t) %in% c("0t_bs2", "24t_bs2", "48t_bs2", "54t_bs2", "60t_bs2", "66t_bs2", "72t_bs2", "96t_bs2", "120t_bs2", "144t_bs2", "168t_bs2"))
tpoints_c_2 = which(rownames(regen_c) %in% c("0c_bs2", "24c_bs2", "48c_bs2", "54c_bs2", "60c_bs2", "66c_bs2", "72c_bs2", "96c_bs2", "120c_bs2", "144c_bs2", "168c_bs2"))

regen_t_1 = regen_t[tpoints_t_1, ]
regen_c_1 = regen_c[tpoints_c_1, ]
regen_t_1 = as.matrix(regen_t_1)
regen_c_1 = as.matrix(regen_c_1)
print("dim(regen_t_1)")
print(rownames(regen_t_1))
print("dim(regen_c_1)")
print(rownames(regen_c_1))

regen_t_2 = regen_t[tpoints_t_2, ]
regen_c_2 = regen_c[tpoints_c_2, ]
regen_t_2 = as.matrix(regen_t_2)
regen_c_2 = as.matrix(regen_c_2)
print("dim(regen_t_2)")
print(rownames(regen_t_2))
print("dim(regen_c_2)")
print(rownames(regen_c_2))

#regen_t_1 = scale(regen_t_1, center=TRUE, scale = FALSE)
#regen_c_1 = scale(regen_c_1, center=TRUE, scale = FALSE)
#regen_t_2 = scale(regen_t_2, center=TRUE, scale = FALSE)
#regen_c_2 = scale(regen_c_2, center=TRUE, scale = FALSE)

top = 1000
var_t = apply(regen_t_1, 2, var)
var_t = sort(var_t, decreasing=TRUE, index.return = TRUE)
#print(paste("var_t$x", var_t$x[1:10],sep=" "))
#print(paste("var_t$ix", var_t$ix[1:10],sep=" "))
regen_t_1 = regen_t_1[ ,var_t$ix[1:top]]
regen_c_1 = regen_c_1[ ,var_t$ix[1:top]]
regen_t_2 = regen_t_2[ ,var_t$ix[1:top]]
regen_c_2 = regen_c_2[ ,var_t$ix[1:top]]

dim = dim(regen_t_1)
n = dim[1]
p = dim[2]

cov_t_1 = (t(regen_t_1) %*% regen_t_1) / (n-1) 
cov_c_1 = (t(regen_c_1) %*% regen_c_1) / (n-1) 
cov_t_2 = (t(regen_t_2) %*% regen_t_2) / (n-1)
cov_c_2 = (t(regen_c_2) %*% regen_c_2) / (n-1)

cov_t_1 = cov_t_1 / max(cov_t_1)
cov_c_1 = cov_c_1 / max(cov_c_1)
cov_t_2 = cov_t_2 / max(cov_t_2)
cov_c_2 = cov_c_2 / max(cov_c_2)
print(cov_t_1[1:10, 1:10])

myHclustDist = function(x) { hclust(d=x, method="average") }

type = "uu"
pdf(paste("../data/cov/",type,"_t_bs1.pdf", sep=""))
hv <- heatmap.2(cov_t_1,trace="none",scale="none",symm = TRUE,col=greenred(75), hclust = myHclustDist, main = "Utricle treated samples, replicate 1") 
dev.off()
pdf(paste("../data/cov/",type,"_c_bs1.pdf", sep=""))
hv <- heatmap.2(cov_c_1,trace="none",scale="none",col=greenred(75), hclust = myHclustDist, main = "Utricle control samples, replicate 1") 
dev.off()
pdf(paste("../data/cov/",type,"_t_bs2.pdf", sep=""))
hv <- heatmap.2(cov_t_2,trace="none",scale="none",col=greenred(75), hclust = myHclustDist, main = "Utricle treated samples, replicate 2") 
dev.off()
pdf(paste("../data/cov/",type,"_c_bs2.pdf", sep=""))
hv <- heatmap.2(cov_c_2,trace="none",scale="none",col=greenred(75), hclust = myHclustDist, main = "Utricle control samples, replicate 2") 
dev.off()

