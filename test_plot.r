roc_fp = matrix(0, 5, 30)
roc_tp = matrix(0, 5, 30)

outdir = "../data/spca7_sim/"

roc_fp[1, ] = c(0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000 , 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.001028807, 0.003086420, 0.007201646, 0.010288066, 0.011316872, 0.011316872)

 roc_tp[1, ] = c(0.0000000, 0.4333333, 0.6333333, 0.7666667, 0.9000000, 0.9000000, 0.9333333, 0.9666667, 
 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000)

roc_fp[2, ] = c(0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.003086420, 0.005144033, 0.007201646, 0.008230453, 0.010288066, 0.013374486, 0.014403292, 0.014403292, 0.015432099, 0.017489712)
 
 roc_tp[2, ] = c(0.0000000, 0.2666667, 0.5333333, 0.6333333, 0.7000000, 0.8000000, 0.8666667, 0.8666667, 0.8666667, 0.9000000, 0.9000000, 0.9000000, 0.9000000, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 1.0000000, 1.0000000, 1.0000000, 1.0000000)
 
 roc_fp[3, ] = c(0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.001028807, 0.003086420, 0.003086420, 0.005144033, 0.007201646, 0.007201646, 0.007201646, 0.006172840, 0.010288066, 0.012345679, 0.013374486, 0.014403292, 0.013374486)

 roc_tp[3, ] = c(0.0000000, 0.3333333, 0.6666667, 0.7666667, 0.8333333, 0.8666667, 0.9000000, 0.9000000, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 0.9666667, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000)

roc_fp[4,] = c(0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.003086420, 0.004115226, 0.004115226, 0.006172840, 0.006172840, 0.006172840)

roc_tp[4, ] = c(0.0000000, 0.2666667, 0.5666667, 0.6666667, 0.8333333, 0.8333333, 0.8333333, 0.8666667, 0.9000000, 0.9000000, 0.9000000, 0.9666667, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000)

roc_fp[5, ] = c(0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.001028807, 0.001028807, 0.001028807, 0.002057613, 0.002057613, 0.002057613, 0.002057613, 0.002057613, 0.003086420, 0.003086420, 0.003086420, 0.003086420, 0.004115226)

roc_tp[5, ] = c(0.0000000, 0.3666667, 0.5333333, 0.6333333, 0.6333333, 0.6666667, 0.8000000, 0.8333333, 0.8333333, 0.8333333, 0.8666667, 0.8666667, 0.9333333, 0.9666667, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000)

#Plot ROC curves for each PC
for (j in 1:1) { 
	#pdf(paste(outdir, "sim2_roc_pc", j,".pdf",sep=""))
	#xaxt="n"
	plot((roc_fp[j, ]), (roc_tp[j, ]), cex=1, pch=16, type="l", xlim=c(0, 0.00001), ylim=c(0,1), xlab ="False positive rate", ylab="True positive rate", main="ROC curves for 2 PCs", cex.lab=1, cex.axis=1)
	#dev.off() 
}