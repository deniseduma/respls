rm(list=ls())
setwd("/home/duma/HairCell/haircell")

t1 = readLines("../data/nscumcomp_t/nscumcomp_genes.txt.diff.1")
t2 = readLines("../data/nscumcomp_t/nscumcomp_genes.txt.diff.2")
inter = intersect(t1, t2)

print("length(t1)")
print(length(t1))
print("length(t2)")
print(length(t2))
print("length(inter)")
print(length(inter))

writeLines(inter, "../data/nscumcomp_t/nscumcomp_genes.txt.diff.12")
