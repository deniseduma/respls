PlotReproduce <- function(fn, v1, v2, n1, n2){
	tiff(file=fn) 
	
	plot(log2(v1), log2(v2), xlab=paste("log2(", n1, ")", sep=""), ylab=paste("log2(", n2, ")", sep=""), pch=19, font=2)
	#v1[v1==0] <- 0.00001
	#v2[v2==0] <- 0.00001
	#v1 <- log2(v1)
	#v2 <- log2(v2)
	r = cor(v1, v2)
	text(2, round(max(log2(v1), log2(v2))*0.9), paste("r = ", round(r, 2), sep=""))
	
	dev.off()
}

GetDESeq_DEG <- function(counts.table, conds, cd1, cd2, fn_normalized, fn_disp, fn_res, fn_res.sorted){
	cds <- newCountDataSet(counts.table, conds)
	cds <- estimateSizeFactors(cds)
	sizeFactors(cds)

	# Get normalized counts
	normTable <- counts(cds, normalized=TRUE)
	# Alternatively, could do the following:
	#{
	#	sf <- sizeFactors(cds)
	#	normTable <- sapply(1:ncol(counts.table), function(i) counts.table[,i]/sf[i])
	#}
	write.table(normTable, file=fn_normalized, sep="\t", quote = FALSE)

	# estimate dispersion (in order to estimate variance for the DE test later)
	cds <- estimateDispersions(cds)
	str(fitInfo(cds))
	tiff(file=fn_disp)
	plotDispEsts(cds)
	dev.off()

	# perform test
	res <- nbinomTest(cds, cd1, cd2)
	save(res, file=fn_res)

	res.high <- res[res$baseMean>0,]
	res.sorted <- res.high[order(res.high$pval),]
	head(res.sorted)
	write.table(res.sorted, file=fn_res.sorted, row.names=F, sep="\t", quote = FALSE)
}

GetDESeq_DEG.blind <- function(counts.table, conds, cd1, cd2, fn_normalized, fn_disp, fn_res, fn_res.sorted){
	cds <- newCountDataSet(counts.table, conds)
	cds <- estimateSizeFactors(cds)
	sizeFactors(cds)

	# Get normalized counts
	normTable <- counts(cds, normalized=TRUE)
	# Alternatively, could do the following:
	#{
	#	sf <- sizeFactors(cds)
	#	normTable <- sapply(1:ncol(counts.table), function(i) counts.table[,i]/sf[i])
	#}
	write.table(normTable, file=fn_normalized, sep="\t", quote = FALSE)

	# estimate dispersion (in order to estimate variance for the DE test later)
	cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only")
	str(fitInfo(cds))
	tiff(file=fn_disp)
	plotDispEsts(cds)
	dev.off()

	# perform test
	res <- nbinomTest(cds, cd1, cd2)
	save(res, file=fn_res)

	res.high <- res[res$baseMean>0,]
	res.sorted <- res.high[order(res.high$pval),]
	head(res.sorted)
	write.table(res.sorted, file=fn_res.sorted, row.names=F, sep="\t", quote = FALSE)
}




DESeq.variance <- function(miu, s, alpha){
	return(s*miu + alpha*(s*s)*(miu*miu))
}

# --------------------------
library(DESeq)

# --------------------------
# Read in counts
setwd("/home/duma/HairCell/haircell")
#ki_M1 <- read.table("../tophat_atxn1ki/ki_M/htseq_count_ensgene.txt", sep="\t", header=F)
#wt_M1 <- read.table("../tophat_atxn1ki/wt_M/htseq_count_ensgene.txt", sep="\t", header=F)
#ki_M2 <- read.table("s8_v2/tophat/ens.nsorted.count.txt", sep="\t", header=F)
#wt_M2 <- read.table("s3_v2/tophat/ens.nsorted.count.txt", sep="\t", header=F)
cc_regen <- read.table("../data/cc_regeneration.csv", sep=",", header=T)

# combine all counts
#ki <- merge(ki_M1, ki_M2, by="V1")
#wt <- merge(wt_M1, wt_M2, by="V1")
#c1 <- merge(ki, wt, by="V1")
#colnames(c1) <- c("gene_ID", "ki_M1", "ki_M2", "wt_M1", "wt_M2")
c1 = cc_regen
print("str(c1)")
print(str(c1))

# Remove extra rows that's not in the feature
#extras <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
#c1 <- subset(c1, !c1[,1] %in% extras)

# Prepare raw counts
counts.table <- c1[,-1]
print("dim(counts.table)")
print(dim(counts.table))
rownames(counts.table) <- c1[,1]
#print("colnames(counts.table)")
#print(colnames(counts.table))
#write.table(counts.table, file="DESeq/counts.table.txt", sep="\t", quote=F)

# plots reproducability graph
#PlotReproduce("DESeq/CorrGraph_kiM1vsKiM2.tiff", counts.table$ki_M1, counts.table$ki_M2, "ki_M1", "ki_M2")
#PlotReproduce("DESeq/CorrGraph_wtM1vswtM2.tiff", counts.table$wt_M1, counts.table$wt_M2, "wt_M1", "wt_M2")
#PlotReproduce("DESeq/CorrGraph_kiM1vswtM1.tiff", counts.table$ki_M1, counts.table$wt_M1, "ki_M1", "wt_M1")
#PlotReproduce("DESeq/CorrGraph_kiM2vswtM2.tiff", counts.table$ki_M2, counts.table$wt_M2, "ki_M2", "wt_M2")

# plots reproducability graph on normalized count
#norm.counts <- read.table("DESeq/norm_KIvsWT_HTseq_count_normalized.txt", sep="\t", header=T)
#PlotReproduce("DESeq/norm_CorrGraph_kiM1vsKiM2.tiff", norm.counts$ki_M1, norm.counts$ki_M2, "norm.ki_M1", "norm.ki_M2")
#PlotReproduce("DESeq/norm_CorrGraph_wtM1vswtM2.tiff", norm.counts$wt_M1, norm.counts$wt_M2, "norm.wt_M1", "norm.wt_M2")
#PlotReproduce("DESeq/norm_CorrGraph_kiM1vswtM1.tiff", norm.counts$ki_M1, norm.counts$wt_M1, "norm.ki_M1", "norm.wt_M1")
#PlotReproduce("DESeq/norm_CorrGraph_kiM2vswtM2.tiff", norm.counts$ki_M2, norm.counts$wt_M2, "norm.ki_M2", "norm.wt_M2")
PlotReproduce("../DESeq/cc_0c_vs_cc_24c.pdf", counts.table$cc_0c, counts.table$cc_24c, "cc_0c", "cc_24c")
PlotReproduce("../DESeq/cc_72c_vs_cc_96c.pdf", counts.table$cc_72c, counts.table$cc_96c, "cc_72c", "cc_96c")
PlotReproduce("../DESeq/cc_144c_vs_cc_168c.pdf", counts.table$cc_144c, counts.table$cc_168c, "cc_144c", "cc_168c")

PlotReproduce("../DESeq/cc_0t_vs_cc_24t.pdf", counts.table$cc_0t, counts.table$cc_24t, "cc_0t", "cc_24t")
PlotReproduce("../DESeq/cc_72t_vs_cc_96t.pdf", counts.table$cc_72t, counts.table$cc_96t, "cc_72t", "cc_96t")
PlotReproduce("../DESeq/cc_144t_vs_cc_168t.pdf", counts.table$cc_144t, counts.table$cc_168t, "cc_144c", "cc_168c")


# ------------------------------------------
#	Plot heatmap
# ------------------------------------------
#head(counts.table)
#conds <- as.factor(c("KI", "KI", "WT", "WT"))
conds <- as.factor(c("cc_0c","cc_0c_bs2","cc_0t","cc_0t_bs2","cc_24c","cc_24c_bs2","cc_24t","cc_24t_bs2","cc_48c","cc_48c_bs2","cc_48t","cc_48t_bs2","cc_72c","cc_72c_bs2","cc_72t","cc_72t_bs2","cc_96c","cc_96c_bs2","cc_96t","cc_96t_bs2","cc_120c","cc_120c_bs2","cc_120t","cc_120t_bs2","cc_144c","cc_144c_bs2","cc_144t","cc_144t_bs2","cc_168c","cc_168c_bs2","cc_168t","cc_168t_bs2"))
cds <- newCountDataSet(counts.table, conds)
cds <- estimateSizeFactors(cds)
print("sizeFactors(cds)")
print(sizeFactors(cds))
cds <- estimateDispersions(cds)

# To avoild -Inf fold changes caused in genes with very low coun, variance stabilizing helps to moderate the fold change estimate - so could plot or cluster
vsd <- getVarianceStabilizedData(cds)
#write.table(vsd, file="DESeq/vsd.txt", sep="\t", quote=F)

dists <- dist(t(vsd))
#tiff(file="DESeq/sample.heatmap2.tif")
pdf(file="../DESeq/heatmap2.pdf")
heatmap( as.matrix(dists), symm=TRUE, scale="none", margins=c(10,10),
	col = colorRampPalette(c("darkblue","white"))(100),
	#labRow = paste( pData(cds)$condition, pData(cds)$libType ) )
	labRow = conditions(cds) )
dev.off()
stop("That's it.")

#  On normalized data
normTable <- counts(cds, normalized=TRUE)
normTable <- read.table("norm_KIvsWT_HTseq_count_normalized.txt", sep="\t", header=T)
head(normTable)
cU <- cor(normTable)
hM <- format(cU, 2)
heatmap.2(cU, Rowv=FALSE, symm=TRUE, col=rev(heat.colors(16)), 
             distfun=function(c) as.dist(1 - c), trace="none", 
             cellnote=hM)

heatmap.2(cU, symm=TRUE, col=rev(heat.colors(16)), 
             distfun=function(c) as.dist(1 - c), trace="none", 
             cellnote=hM)


# Plot heatmap on the significantly differentiate genes
sig.dat <- read.table("DESeq/DESeq.SigGenes.normTable.txt", sep="\t", header=T)
sig.vsd <- subset(vsd, row.names(vsd) %in% sig.dat[,1])
tiff(file="DESeq/sig.heatmap.tif")
heatmap(sig.vsd, col = colorRampPalette(c("darkblue","white"))(100), scale="none", labRow=sig.dat[,2])
dev.off()


# ------------------------------------------
#	DE Analysis
# ------------------------------------------
# Make the data set required for DE test
conds <- as.factor(c("KI", "KI", "WT", "WT"))

cd1 = "KI"
cd2 = "WT"
fn_normalized = "DESeq/norm_KIvsWT_HTseq_count_normalized.txt"
fn_disp =  "DESeq/norm_dispersionplot.tif"
fn_res = "DESeq/norm_KIvsWT_res.saved"
fn_res.sorted ="DESeq/norm_KIvsWT_DEseq.txt"

GetDESeq_DEG(counts.table, conds, cd1, cd2, fn_normalized, fn_disp, fn_res, fn_res.sorted)





# --------------
# Hyojin's code
library(DESeq)
setwd("/home/hjkang/Project/NGS/RNASeq/Atxn1KO/deseq/")
#d = read.table(file="htseq_count_refgene_mm9.txt", sep="\t", header=T)
countsTable = cbind(d$WT1, d$WT2, d$KO1, d$KO2)
dim(countsTable)
colnames(countsTable) <- c("WT1", "WT2", "KO1", "KO2")
rownames(countsTable) <- d$X
head(countsTable)
conds <- as.factor(c("WT","WT", "KO","KO"))
cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
normTable <- counts(cds, normalized=TRUE)
write.table(normTable, file="htseq_count_normalized.txt", sep="\t", quote = FALSE)
cds <- estimateDispersions(cds)
#perform test
res <- nbinomTest(cds, "WT", "KO")
res.high <- res[res$baseMean>0,]
res.sorted <- res.high[order(res.high$pval),]
head(res.sorted)
#write.table(res.sorted, file="Atxn1KO_deseq_ensg.txt", row.names=F, sep="\t", quote = FALSE)
