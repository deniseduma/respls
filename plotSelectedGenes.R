rm(list=ls())
setwd("/home/duma/HairCell/haircell")

mytitle = "Selected genes"

regen_cc = read.table(paste("../data/cc_regeneration.csv", sep = ""), sep=",", header=T, row.names=1)
apoints_cc = readLines(paste("../data/cc_timepoints.txt", sep=""))
tpoints_cc = readLines(paste("../data/cc_timepoints_t.txt", sep=""))

regen_uu = read.table(paste("../data/uu_regeneration.csv", sep = ""), sep=",", header=T, row.names=1)
apoints_uu = readLines(paste("../data/uu_timepoints.txt", sep=""))
tpoints_uu = readLines(paste("../data/uu_timepoints_t.txt", sep=""))

colnames(regen_uu) = apoints_uu

diff = readLines("../data/nscumcomp_t/nscumcomp_genes_diff.txt")
print("length(diff)")
print(length(diff))

regen_cc = t(regen_cc)
regen_uu = t(regen_uu)

idx_cc = which(colnames(regen_cc) %in% diff)
idx_uu = which(colnames(regen_uu) %in% diff)

timepts_cc = which(rownames(regen_cc) %in% tpoints_cc)
timepts_uu = which(rownames(regen_uu) %in% tpoints_uu)

pdf("../data/nscumcomp_t/nscumcomp_gene_diff.pdf")
par(mfrow=c(1, 2))
matplot(timepts_cc, regen_cc[timepts_cc, idx_cc], type="l", pch=16, main = mytitle, xlab = "Cochlea timepts", ylab = "Gene expression")
matplot(timepts_uu, regen_uu[timepts_uu, idx_uu], type="l", pch=16, main = mytitle, xlab = "Utricle timepts", ylab = "Gene expression")
dev.off()

print("regen_cc[timepts_cc, idx_cc]")
print(regen_cc[timepts_cc, idx_cc])
print("regen_uu[timepts_uu, idx_uu]")
print(regen_uu[timepts_uu, idx_uu])

