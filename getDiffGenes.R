#!/usr/bin/Rscript

rm(list=ls())
setwd("/home/duma/HairCell/haircell")

args=commandArgs(TRUE)
version = args[1] 
threshold = args[2]

#t = readLines("../data/nscumcomp_t/nscumcomp_genes.txt")
#c = readLines("../data/nscumcomp_c/nscumcomp_genes.txt")
t = readLines(paste("../data/nscumcomp_t/", threshold, "TimesGenes.txt.t.", version, sep=""))
c = readLines(paste("../data/nscumcomp_c/", threshold, "TimesGenes.txt.c.", version, sep=""))
diff = setdiff(t, c)

print("length(t)")
print(length(t))
print("length(c)")
print(length(c))
print("length(diff)")
print(length(diff))

#writeLines(diff, "../data/nscumcomp_t/nscumcomp_genes.txt.diff")
writeLines(diff, paste("../data/nscumcomp_t/", threshold, "TimesGenes.txt.diff.", version, sep=""))
