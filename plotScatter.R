rm(list=ls())

U = read.table('../data/U_Multi_CT.csv', header=T, sep=",", row.names=1, stringsAsFactors=F)
U = as.matrix(U)

V1 = read.table('../data/V1_Multi_CT.csv', header=T, sep=",", row.names=1, stringsAsFactors=F)
V1 = as.matrix(V1)

V2 = read.table('../data/V2_Multi_CT.csv', header=T, sep=",", row.names=1, stringsAsFactors=F)
V2 = as.matrix(V2)

X1 = read.table('../data/X1_CT.csv', header=T, sep=",", row.names=1, stringsAsFactors=F)
X1 = as.matrix(X1)

Y0 = read.table('../data/Y0_CT.csv', header=T, sep=",", row.names=1, stringsAsFactors=F)
Y0 = as.matrix(Y0)

#pcs=myPCA(train=ludata[index,],input=ludata[index,],label=label,i=1,j=2,shape=16,lsize=lsize)

pdf("../figs/X1_Cols_CT.pdf")
boxplot(as.data.frame(X1), ylab="Columns of X1", main="X1")
dev.off()

pdf("../figs/Y0_Cols_CT.pdf")
boxplot(as.data.frame(Y0), ylab="Columns of Y0", main="Y0")
dev.off()

pdf("../figs/U_Multi_Cols_CT.pdf")
boxplot(as.data.frame(U), ylab="Columns of U", main="U_Multi")
dev.off()

#barplots for the columns of U
for(i in 1:(ncol(U))) {
	pdf(paste("../figs/U_Multi-PC", i, "_CT.pdf", sep=""))
	barplot(U[,i], main=paste("U_Multi-PC", i, sep=""))
	dev.off()
}

#colors = rep(c("red", "blue", "pink", "green", "cyan", "magenta", "grey", "orange"), 4);
colors = c(rep("green", 8), rep("red", 8), rep("green", 8), rep("red", 8));

lsize = rep(c(0,24,48,72,96,120,144,168), 4)
lsize = lsize/max(lsize)*6 + 2

print(nrow(V1))
print(ncol(V1))

#scatter plots of components of V1
for(i in 1:(ncol(V1)-1)){
  for(j in (i+1):ncol(V1)){
    pdf(paste("../figs/V1_Multi-PC", i, "-", j, "_CT.pdf", sep=""))
    plot(V1[,i], V1[,j], pch=19, cex = lsize, col = colors, xlab=paste("PC", i), ylab=paste("PC", j), main="Cochlea") 
    legend(2000, 9.5, c("Control", "Treated"), lty=c(1,1), lwd=c(2.5,2.5), col=c("green", "red"))
    #abline(0,1,col="red")
    dev.off()
  }
}

ct_bs3 = c(24, 72, 96, 120, 168)
ct_bs2 = c(0,24,48,54,60,66,72,96,120,144,168) 
#lsize = c(rep(ct_bs2, 2), c(0, 0), rep(ct_bs2, 2), rep(ct_bs3, 2), c(24, 24))
lsize = c(rep(ct_bs2, 2), rep(ct_bs2, 2))
lsize = lsize/max(lsize)*6+1

#ct_bs3 = c("blue", "green", "cyan", "magenta", "orange")
#ct_bs2 = c("red", "blue", "pink", "darkgoldenrod4", "brown", "blueviolet", "green", "cyan", "magenta", "grey", "orange");
#colors = c(rep(ct_bs2, 2), c("red", "red"), rep(ct_bs2, 2), rep(ct_bs3, 2), c("blue", "blue"))
#colors = c(rep("green", 11), rep("red", 11), "green", "red", rep("green", 11), rep("red", 11), rep("green", 5), rep("red", 5), "green", "red")
colors = c(rep("green", 11), rep("red", 11), rep("green", 11), rep("red", 11))

print(nrow(V2))
print(ncol(V2))

#scatter plots of components of V2
for(i in 1:(ncol(V2)-1)){
  for(j in (i+1):ncol(V2)){
    pdf(paste("../figs/V2_Multi-PC", i, "-", j, "_CT.pdf", sep=""))
    plot(V2[,i], V2[,j], pch=19, cex = lsize, col = colors, xlab=paste("PC", i), ylab=paste("PC", j), main="Utricle")
    legend(2000, 9.5, c("Control", "Treated"), lty=c(1,1), lwd=c(2.5,2.5), col=c("green", "red"))
    #abline(0,1,col="red")
    dev.off()
  }
}


