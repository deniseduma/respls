
printHeatmap <- function(X) {

	#Color scheme for heatmap
	my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
	#(Optional) defines the color breaks manually for a "skewed" color transition
	#col_breaks = c(seq(-1,0,length=100),  # for red
	#  seq(0,0.8,length=100),              # for yellow
	#  seq(0.8,1,length=100))              # for green

	#Heatmap of input matrix
	#pdf(paste(outdir, "x.pdf", sep=""))
	png(paste(outdir,"x.png",sep=""),    # create PNG for the heat map        
	width = 5*300,        # 5 x 300 pixels
	height = 5*300,
	res = 300,            # 300 pixels per inch
	pointsize = 8)
	
	myHclustDist = function(x) { hclust(d=x, method="average") }
	hv<-heatmap.2(X,trace="none",scale="none",margins=c(12,9),col=my_palette,dendrogram="col",hclust = myHclustDist,main = "X") 
	
	dev.off()

	#Heatmap of input matrix cov
	#pdf(paste(outdir, "xx.pdf", sep=""))
	#png(paste(outdir,"xx.png",sep=""),    # create PNG for the heat map        
	#width = 5*300,        # 5 x 300 pixels
	#height = 5*300,
	#res = 300,            # 300 pixels per inch
	#pointsize = 8)
	
	#myHclustDist = function(x) { hclust(d=x, method="average") }
	#hv<-heatmap.2(t(regen)%*%regen,trace="none",scale="none",margins=c(12,9),col=my_palette,hclust = myHclustDist,main = "XX") 
	
	#dev.off()

}
