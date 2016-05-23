rm(list=ls())
setwd("/home/duma/HairCell/haircell")

regen <- read.table("../data/uu_regeneration.csv", sep=",", header=T, row.names=1)
apoints = readLines("../data/uu_timepoints.txt")
cpoints = readLines("../data/uu_timepoints_c.txt")
tpoints = readLines("../data/uu_timepoints_t.txt")

colnames(regen) = apoints
GENES = rownames(regen)

print("dim(regen)")
print(dim(regen))
print("length(rownames(regen))")
print(length(rownames(regen)))
print("length(colnames(regen)")
print(length(colnames(regen)))
print("head(colnames(regen))")
print(head(colnames(regen)))

control = regen[, cpoints]
print("dim(control)")
print(dim(control))

treated = regen[, tpoints]
print("dim(treated)")
print(dim(treated))

#only keep genes whose max foldchange  >= 0.5
rowMax = apply(as.matrix(treated), 1, max)
print("head(rowMax)")
print(head(rowMax))
regen0 = regen[rowMax >= 0.5, ]
print("dim(regen0)")
print(dim(regen0))
#print("head(regen1)")
#print(head(regen1))

GENES0 = rownames(regen0)
print("length(GENES0)")
print(length(GENES0))

myFilter = function(gene, tc) {
	test = kruskal.test(gene ~ tc)
	return(test$p.value)
}

tc = colnames(regen0)
tc[which(tc %in% cpoints)] = rep("C", 30)
tc[which(tc %in% tpoints)] = rep("T", 30)
tc = factor(tc)
#filter for genes which vary between treated and control
p.values = apply(as.matrix(regen0), 1, function(l) myFilter(l, tc))
regen1 = regen0[p.values <= 0.05, ]
print("dim(regen1)")
print(dim(regen1))

GENES1 = rownames(regen1)
print("length(GENES1)")
print(length(GENES1))

writeLines(GENES1, "../data/uu_by_tc.txt")

tc = as.factor(c(0,0,0,24,24,24,24,48,48,54,54,60,60,66,66,72,72,72,96,96,96,120,120,120,120,144,144,168,168,168))
#filter for genes which vary across time points
p.values = apply(as.matrix(regen1), 1, function(l) myFilter(gene=l[which(names(l) %in% tpoints)], tc) )
regen2 = regen1[p.values <= 0.05, ]
print("dim(regen2)")
print(dim(regen2))

GENES2 = rownames(regen2)
print("length(GENES2)")
print(length(GENES2))

writeLines(GENES2, "../data/uu_by_tc_time.txt")



