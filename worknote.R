### read in the ku file. 
rm(list=ls())
udata = read.table('../data/Ku.txt',header=T,sep="\t",row.names=1,stringsAsFactors=F)
udata = as.matrix(udata)

cdata = read.table('../data/CC_regeneration_all.txt',header=T,sep="\t",row.names=1,stringsAsFactors=F,skip=)
cdata = as.matrix(cdata)

snames=intersect(rownames(udata),rownames(cdata))
udata=udata[snames,]
cdata=cdata[snames,]

ludata = log2(1+udata)
lcdata = log2(1+cdata)


#### PCA on data. 

for(i in 1:(ncol(ludata)-1)){
  for(j in (i+1):ncol(ludata)){
    pdf(paste("../output/lu",i,"-",j,".pdf",sep=""))
    plot(ludata[,i],ludata[,j],pch=16)
    abline(0,1,col="red")
    dev.off()
  }
}

for(i in 1:(ncol(lcdata)-1)){
  for(j in (i+1):ncol(lcdata)){
    pdf(paste("../output/lc",i,"-",j,".pdf",sep=""))
    plot(lcdata[,i],lcdata[,j],pch=16)
    abline(0,1,col="red")
    dev.off()
  }
}

pdf("../output/boxplot-cochlea.pdf")
boxplot(as.data.frame(lcdata),ylab="log(1+X)",main="cochlea")
dev.off()
pdf("../output/boxplot-utricle.pdf")
boxplot(as.data.frame(ludata),ylab="log(1+X)",main="utricle")
dev.off()

############   

ss= apply(ludata,1,sd)
index = names(sort(ss,decreasing=T)[1:1000]);


label = c(0,0,0,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,1)
label[label==0]="red"
label[label==1]="green"

lsize = rep(c(0,24,48,54,60,66,72,96,120,144,168),c(6,8,4,4,4,4,6,6,8,4,6))
lsize = lsize/max(lsize)*5+1

pdf("../output/Utricle-PCA1-2-t1000.pdf")
pcs=myPCA(train=ludata[index,],input=ludata[index,],label=label,i=1,j=2,shape=16,lsize=lsize)
dev.off()

pdf("../output/Utricle-PCA1-3-t1000.pdf")
myPCA(train=ludata[index,],input=ludata[index,],label=label,i=1,j=3,shape=16,lsize=lsize)
dev.off()


pdf("../output/Utricle-PCA1-4-t1000.pdf")
myPCA(train=ludata[index,],input=ludata[index,],label=label,i=1,j=4,shape=16,lsize=lsize)
dev.off()


pdf("../output/Utricle-PCA2-3-t1000.pdf")
myPCA(train=ludata[index,],input=ludata[index,],label=label,i=2,j=3,shape=16,lsize=lsize)
dev.off()

pdf("../output/Utricle-PCA2-4-t1000.pdf")
myPCA(train=ludata[index,],input=ludata[index,],label=label,i=2,j=4,shape=16,lsize=lsize)
dev.off()

pdf("../output/Utricle-PCA3-4-t1000.pdf")
myPCA(train=ludata[index,],input=ludata[index,],label=label,i=3,j=4,shape=16,lsize=lsize)
dev.off()


pdf("../output/Utricle-PC1.pdf")
barplot(pcs$rotation[,1],main="Utricle-PC1")
dev.off()

pdf("../output/Utricle-PC2.pdf")
barplot(pcs$rotation[,2],main="Utricle-PC2")
dev.off()

pdf("../output/Utricle-PC3.pdf")
barplot(pcs$rotation[,3],main="Utricle-PC3")
dev.off()


write(index[abs(pcs$rotation[,1])>0.05], file="../output/Utri-PC1.txt")
write(index[abs(pcs$rotation[,2])>0.05], file="../output/Utri-PC2.txt")
write(index[abs(pcs$rotation[,3])>0.05], file="../output/Utri-PC3.txt")


#######################  PCA on Cochlea

ss= apply(lcdata,1,sd)
index = names(sort(ss,decreasing=T)[1:1000]);

label = c(0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1)
label[label==0]="red"
label[label==1]="green"

lsize = rep(c(0,24,48,72,96,120,144,168),4)
lsize = lsize/max(lsize)*5+1

pdf("../output/Cochlea-PCA1-2-t1000.pdf")
pcs=myPCA(train=lcdata[index,],input=lcdata[index,],label=label,i=1,j=2,shape=16,lsize=lsize)
dev.off()

pdf("../output/Cochlea-PCA1-3-t1000.pdf")
myPCA(train=lcdata[index,],input=lcdata[index,],label=label,i=1,j=3,shape=16,lsize=lsize)
dev.off()


pdf("../output/Cochlea-PCA1-4-t1000.pdf")
myPCA(train=lcdata[index,],input=lcdata[index,],label=label,i=1,j=4,shape=16,lsize=lsize)
dev.off()


pdf("../output/Cochlea-PCA2-3-t1000.pdf")
myPCA(train=lcdata[index,],input=lcdata[index,],label=label,i=2,j=3,shape=16,lsize=lsize)
dev.off()

pdf("../output/Cochlea-PCA2-4-t1000.pdf")
myPCA(train=lcdata[index,],input=lcdata[index,],label=label,i=2,j=4,shape=16,lsize=lsize)
dev.off()

pdf("../output/Cochlea-PCA3-4-t1000.pdf")
myPCA(train=lcdata[index,],input=lcdata[index,],label=label,i=3,j=4,shape=16,lsize=lsize)
dev.off()

pdf("../output/Cochlea-PC1.pdf")
barplot(pcs$rotation[,1],main="Cochlea-PC1")
dev.off()

pdf("../output/Cochlea-PC2.pdf")
barplot(pcs$rotation[,2],main="Cochlea-PC2")
dev.off()

pdf("../output/Cochlea-PC3.pdf")
barplot(pcs$rotation[,3],main="Cochlea-PC3")
dev.off()

write(index[abs(pcs$rotation[,1])>0.05], file="../output/Coch-PC1.txt")
write(index[abs(pcs$rotation[,2])>0.05], file="../output/Coch-PC2.txt")
write(index[abs(pcs$rotation[,3])>0.05], file="../output/Coch-PC3.txt")








