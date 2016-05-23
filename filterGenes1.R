rm(list=ls())
setwd("/home/duma/HairCell/haircell")

type = "uu"

## read cc
regen = read.table(paste("../data/", type, "_regen.csv", sep = ""), sep=",", header=T, row.names=1)
apoints = readLines(paste("../data/", type, "_timepoints.txt", sep=""))
cpoints = readLines(paste("../data/", type, "_timepoints_c.txt", sep=""))
tpoints = readLines(paste("../data/", type, "_timepoints_t.txt", sep=""))
colnames(regen) = apoints

regen = t(as.matrix(regen))
regen = scale(regen, center=T, scale=F)
print("dim(regen)")
print(dim(regen))

genes = colnames(regen)

control = regen[cpoints, ]
print("dim(control)")
print(dim(control))

treated = regen[tpoints, ]
print("dim(treated)")
print(dim(treated))

#only keep genes whose max foldchange  >= 0.5
#colMax = apply(treated, 2, max)
#print("head(rowMax)")
#print(head(rowMax))
#regen = regen[, colMax >= 0.5]
#print("dim(regen)")
#print(dim(regen))

lmp = function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f = summary(modelobject)$fstatistic
    p = pf(f[1], f[2], f[3], lower.tail=F)
    attributes(p) = NULL
    return(p)
}

myFilter = function(gene, tc, timepts) {
	test = lm(gene ~ tc:timepts)
	return(lmp(test))
}

## We filter for genes which vary between treated and control samples
tc = rownames(regen)
if (identical(type, "cc")) {
	tc[which(tc %in% cpoints)] = rep("C", 16)
	tc[which(tc %in% tpoints)] = rep("T", 16)
} else if (identical(type, "uu")) {
	tc[which(tc %in% cpoints)] = rep("C", 30)
	tc[which(tc %in% tpoints)] = rep("T", 30)
}
tc = factor(tc)
print("tc")
print(t(tc))

## We also filter for genes which vary across time points
if (identical(type, "cc")) {
	timepts = c(0,0,0,0,24,24,24,24,48,48,48,48,72,72,72,72,96,96,96,96,120,120,120,120,144,144,144,144,168,168,168,168)
} else if (type == "uu") {
	timepts = c(rep(0, 6), rep(24,8), rep(48, 4), rep(54,4), rep(60,4), rep(66,4), rep(72,6), rep(96, 6), rep(120, 8), rep(144, 4), rep(168,6))
}
print("timepts")
print(t(timepts))

## Do the filtering
#p.values = apply(as.matrix(regen1), 1, function(l) myFilter(gene=l[which(names(l) %in% tpoints)], tc, timepts) )
p.values = apply(regen, 2, function(l) myFilter(l, tc, timepts) )
p.values = p.adjust(p.values, method = "BH")
print("length(p.values)")
print(length(p.values))

regen = regen[, which(p.values <= 0.05)]
print("dim(regen)")
print(dim(regen))

genes = colnames(regen)
print("length(genes)")
print(length(genes))

writeLines(genes, paste("../data/inter_by_both1.txt.uu", sep=""))

