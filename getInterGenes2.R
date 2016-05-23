rm(list=ls())
setwd("/home/duma/HairCell/haircell")

threshold = 80
version = 6
oversion = 3

d1 = readLines(paste("../data/nscumcomp_t/", threshold, "TimesGenes.txt.diff.", version, sep=""))
d2 = readLines(paste("../data/nscumcomp_t/", threshold, "TimesGenes.txt.diff.", oversion, sep=""))
inter = intersect(d1, d2)

print("length(d1)")
print(length(d1))
print("length(d2)")
print(length(d2))
print("length(inter)")
print(length(inter))

writeLines(inter, paste("../data/nscumcomp_t/", threshold, "TimesGenes.txt.diff.inter.", oversion, ".", version, sep=""))
