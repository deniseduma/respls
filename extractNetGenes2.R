rm(list=ls())
setwd("/home/duma/HairCell/haircell")

type = "cc"

geneFile = "../data/STRING.txt"
inFile = paste("../data/", type,"_regen.csv", sep = "")
#outFile = paste("../data/", type, "_regen_net.csv", sep = "")
outFile = "../data/inter_by_both2_net.txt"

geneSubset = readLines("../data/inter_by_both2.txt")

X = read.table(inFile, sep=",", header = TRUE, row.names = 1) 
apoints = readLines(paste("../data/", type, "_timepoints.txt", sep=""))
if (identical(type, "uu")) {
	colnames(X) = apoints
}

X = t(X[geneSubset, ])
print("dim(X)")
print(dim(X))

GENES = readLines(geneFile)
GENES = strsplit(GENES, "\t")
print("Number of interactions")
print(length(GENES))
print("Example interactions")
print(head(GENES))

V1 = sapply(GENES, function(l) l[1])
V2 = sapply(GENES, function(l) l[2])
VV = c(V1, V2)
I = unique(VV)
print("Number of unique genes")
print(length(I))

X = X[, which(colnames(X) %in% I)]
print("dim(X)")
print(dim(X))

#write.table(X, file = outFile, quote = FALSE)
writeLines(colnames(X), outFile)
