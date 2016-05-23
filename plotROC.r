
tpr = matrix(c(0.13,0.175,0.283,0.45,0.733,0.875,1,1,1,1),2,5)
print("tpr")
print(tpr)

fpr = matrix(0,2,5)
print("fpr")
print(fpr)

for (j in 1:2) {
	dev.new()
	plot(fpr[j,], tpr[j,], cex=1, pch=16, type="l", xlim=c(0,0.001), ylim=c(0,1), xlab ="fpr", ylab="tpr", main=paste("ROC curve for SPCA PC ",j," n=10,20,50,100,1000; p=100",sep=""), cex.lab=1.25, cex.axis=1.25)
	points(fpr[j,],tpr[j,])
}	


