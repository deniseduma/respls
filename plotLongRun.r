
#Read consensus matrices from file
consensus_t = read.table("../data/cc/consensus_k8_t.txt.lr", quote = "")
consensus_c = read.table("../data/cc/consensus_k8_c.txt.lr", quote = "")
print("dim(consensus_t)")
print(dim(consensus_t))
print("dim(consensus_c)")
print(dim(consensus_c))

#Remove 0 rows and cols
sum1 = apply(consensus_t, 1, sum)
idx1 = which(sum1 == 0)
sum2 = apply(consensus_t, 2, sum)
idx2 = which(sum2 == 0)
consensus_t = consensus_t[-idx1, -idx2]
maxc = apply(consensus_t, 2, max)
print("dim(consensus_t)")
print(dim(consensus_t))
#print("maxc")
#print(maxc)
sum1 = apply(consensus_c, 1, sum)
idx1 = which(sum1 == 0)
sum2 = apply(consensus_c, 2, sum)
idx2 = which(sum2 == 0)
consensus_c = consensus_c[-idx1, -idx2]
maxc = apply(consensus_c, 2, max)
print("dim(consensus_c)")
print(dim(consensus_c))
#print("maxc")
#print(maxc)

#Plot heatmap of consensus matrix
library("gplots")
#myDist = function(x) {dist(x, method="manhattan")}
myDist = function(x) {as.dist(1-x)}
myHclustDist = function(x) { hclust(d=x, method="average") }

#source("unByteCode.R")
#dend = as.dendrogram(myHclustDist(myDist(consensus_t)))
## convert stats:::plotNode from byte-code to interpreted-code
#unByteCodeAssign(stats:::plotNode)
#options("expressions"=5e4)  # increase recursion limit
#plot(dend)

#Treated
pdf("../data/cc/consensus_k8_t.txt.lr.pdf")
heatmap.2(as.matrix(consensus_t), dendrogram="column", symm=TRUE, dist = myDist, hclust = myHclustDist, col=redgreen(75), key=TRUE, density.info="none", trace="none", main=paste("Treated consensus matrix, k=", myncomp,sep=""))
dev.off()
#Control
pdf("../data/cc/consensus_k8_c.txt.lr.pdf")
heatmap.2(as.matrix(consensus_c), dendrogram = "column", symm=TRUE, dist = myDist, hclust = myHclustDist, col=redgreen(75), key=TRUE, density.info="none", trace="none", main=paste("Control consensus matrix, k=", myncomp,sep=""))
dev.off()

