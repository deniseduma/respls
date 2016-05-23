rm(list=ls())
setwd("/home/duma/HairCell/haircell")

inMatrix = "V"
inFile = paste("../data/xpca/", inMatrix, ".csv", sep="")
outFile = paste("../data/xpca/", inMatrix, "_PCs.pdf", sep="")

#Read in the matrix of PCs
X = read.table(inFile, sep=",", header=T, row.names=1)
print("dim(X)") 
print(dim(X)) 

tpoints_cc = readLines(paste("../data/cc_timepoints_t.txt", sep=""))
tpoints_uu = readLines(paste("../data/uu_timepoints_t.txt", sep=""))

label = rownames(X) %in% tpoints_cc
label[label==TRUE]="green"
label[label==FALSE]="red"

lsize = c(rep(0, 2), rep(24,2), rep(48, 2), rep(72, 2), rep(96, 2), rep(120, 2), rep(144,2), rep(168, 2), rep(0, 3), rep(24,4), rep(48, 2), rep(54,2), rep(60,2), rep(66,2), rep(72,3), rep(96, 3), rep(120, 4), rep(144, 2), rep(168,3))
#lsize = lsize/max(lsize)*5+1
lsize = lsize/max(lsize)*1.5+0.5

#scatterplots of the PCs
pdf(outFile)
pairs(X, col=label, pch=16, cex=lsize, main=paste("PCs of ", inMatrix,sep=""))
dev.off()
