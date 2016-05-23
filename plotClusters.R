rm(list=ls())
setwd("/home/duma/HairCell/haircell")

inFile = "../data/xpca/B_K.csv"
outFile = "../data/xpca/B_K.pdf"

#Read in the matrix of PCs
X = read.table(inFile, sep=",", header=T, row.names=1)
print("dim(X)") 
print(dim(X)) 

#scatterplots of the PAs
pdf(outFile)
plot(X, col="red", pch=16, main="Clusters")
dev.off()
