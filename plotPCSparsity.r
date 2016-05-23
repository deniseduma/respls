#!/usr/bin/Rscript
noSteps = 25 

sp = readLines("../data/sanity_check_log/sparsity.txt")
sp = as.numeric(sp)
mymat = matrix(sp[1:6])
for (i in seq(7,150,6)) 
	mymat = cbind(mymat, sp[i:(i+5)])
rownames(mymat) = c("PC1", "PC2", "PC3","PC4","PC5","PC6")
print("mymat")
print(mymat)

pdf("../data/sanity_check_log/sparsity.pdf")
matplot(1:noSteps, t(mymat), type="l", pch=16, col=c(1,2,3,4,5,6), xlab = "l1 penalty rank", ylab = "PC sparsity", main="Sparsity of PCs vs l1 penalty rank") 
legend("topleft", inset=.05, legend=c("PC1", "PC2", "PC3","PC4","PC5","PC6"), pch=16, col=c(1,2,3,4,5,6))
#axis(1, at=1:noSteps)
dev.off()

