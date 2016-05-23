#!/usr/bin/Rscript

library(gplots)

## read cc
regen_cc = read.table(paste("../data/cc_regen.csv", sep = ""), sep=",", header=T, row.names=1)
apoints_cc = readLines(paste("../data/cc_timepoints.txt", sep=""))
cpoints_cc = readLines(paste("../data/cc_timepoints_c.txt", sep=""))
tpoints_cc = readLines(paste("../data/cc_timepoints_t.txt", sep=""))
colnames(regen_cc) = apoints_cc

## read uu
regen_uu = read.table(paste("../data/uu_regen.csv", sep = ""), sep=",", header=T, row.names=1)
apoints_uu = readLines(paste("../data/uu_timepoints.txt", sep=""))
cpoints_uu = readLines(paste("../data/uu_timepoints_c.txt", sep=""))
tpoints_uu = readLines(paste("../data/uu_timepoints_t.txt", sep=""))
colnames(regen_uu) = apoints_uu

regen_cc = as.matrix(regen_cc)
regen_uu = as.matrix(regen_uu)
print("dim(regen_cc)")
print(dim(regen_cc))
print("dim(regen_uu)")
print(dim(regen_uu))

geneSubset = readLines("../data/inter_by_both.txt")
regen_cc = regen_cc[geneSubset, ]
regen_uu = regen_uu[geneSubset, ]

regen_big = cbind(regen_cc, regen_uu)
regen_big = t(regen_big)
print("dim(regen_big)")
print(dim(regen_big))

regen_t =  regen_big[which(rownames(regen_big) %in% c(tpoints_cc, tpoints_uu)), ]
regen_c =  regen_big[which(rownames(regen_big) %in% c(cpoints_cc, cpoints_uu)), ]

genes = colnames(regen_big)

tpoints = c(0, 24, 48, 72, 96, 120, 144, 168)
tpoints_t = which(rownames(regen_t) %in% c("cc_0t", "cc_24t", "cc_48t", "cc_72t", "cc_96t", "cc_120t", "cc_144t", "cc_168t"))
tpoints_c = which(rownames(regen_c) %in% c("cc_0c", "cc_24c", "cc_48c", "cc_72c", "cc_96c", "cc_120c", "cc_144c", "cc_168c"))
#tpoints_t = which(rownames(regen_t) %in% c("cc_0t_bs2", "cc_24t_bs2", "cc_48t_bs2", "cc_72t_bs2", "cc_96t_bs2", "cc_120t_bs2", "cc_144t_bs2", "cc_168t_bs2"))
#tpoints_c = which(rownames(regen_c) %in% c("cc_0c_bs2", "cc_24c_bs2", "cc_48c_bs2", "cc_72c_bs2", "cc_96c_bs2", "cc_120c_bs2", "cc_144c_bs2", "cc_168c_bs2"))

#tpoints_t = which(rownames(regen_t) %in% c("0t", "24t", "48t", "72t", "96t", "120t", "144t", "168t"))
#tpoints_c = which(rownames(regen_c) %in% c("0c", "24c", "48c", "72c", "96c", "120c", "144c", "168c"))
#tpoints_t = which(rownames(regen_t) %in% c("0t_bs2", "24t_bs2", "48t_bs2", "72t_bs2", "96t_bs2", "120t_bs2", "144t_bs2", "168t_bs2"))
#tpoints_c = which(rownames(regen_c) %in% c("0c_bs2", "24c_bs2", "48c_bs2", "72c_bs2", "96c_bs2", "120c_bs2", "144c_bs2", "168c_bs2"))

dimx = dim(regen_t)
n = dimx[1]
p = dimx[2]

regen_t = scale(regen_t, center=TRUE, scale = FALSE)
regen_c = scale(regen_c, center=TRUE, scale = FALSE)

myHclustDist = function(x) { hclust(d=x, method="average") }

iter = 2
for (i in c(1,2,3,4)) {
	
	pdf(paste("../data/cc/v3_noscale_top15_k6_v",iter,"/iter", iter, "_hmap",i,"_log.pdf", sep=""))
	
	#Extract genes to plot
	genesToPlot = readLines(paste("../data/cc/v3_noscale_top15_k6_v",iter,"/cluster", i, ".txt", sep=""))
	print("length(genesToPlot)")
	print(length(genesToPlot))
	
	genesToPlotT = regen_t[tpoints_t, genesToPlot]
	genesToPlotC = regen_c[tpoints_c, genesToPlot]
	genesToPlotT = as.matrix(genesToPlotT)
	genesToPlotC = as.matrix(genesToPlotC)
	foldChange = log2(genesToPlotT / genesToPlotC)
	foldChange[is.infinite(foldChange) | is.nan(foldChange)] = max(foldChange[is.finite(foldChange)], na.rm=TRUE) 
	print("dim(genesToPlotT)")
	#print(dim(genesToPlotT))
	print(foldChange)
	
	hv<-heatmap.2(foldChange,trace="none",scale="none",col=greenred(75), hclust = myHclustDist, main = paste("Fold change, iter ", iter, " cluster ",i, sep="")) 
	
	dev.off()
}


