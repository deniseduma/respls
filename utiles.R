plot.dens<-function(x,y,xlab="",ylab=".Rapp.history",xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),transparency=70,fc=log2(2),cex=1,cex.axis=1){
  
  require(geneplotter)
  
  mdata=cbind(x,y)
  
  
  colMap=heat.colors(100)
  colMap=densCols(mdata, colramp=colorRampPalette(colMap))
  colMap=apply(col2rgb(colMap),2,function(r){rgb(r[1],r[2],r[3], transparency,maxColorValue=255)})
  
  plot(mdata,col=colMap,pch=16,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,cex=cex,cex.axis=cex.axis)
  
  abline(a=fc,b=1)
  abline(a=-fc,b=1)
  
  
} 



myPCA = function(train,input=train,label,shape=10,i,j,lsize=2)
{
  PCAs=prcomp(t(train));
  mapped=t(input)%*%PCAs$rotation;
  plot(mapped[,i],mapped[,j],col=label,pch=shape,cex=lsize,xlab=paste("PC",i),ylab=paste("PC",j))
  #scatterplot3d(mapped[,1:3],color=label,pch=shape)
  text(mapped[,i],mapped[,j],labels=colnames(input),cex=1,offset=0,pos=1)
  print(summary(PCAs))
  return(PCAs)
}
