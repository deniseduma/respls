#!/usr/bin/Rscript

rm(list=ls())
setwd("/home/duma/HairCell/haircell")

#s1 = readLines("../data/inter_by_both1.txt.cc")
#s2 = readLines("../data/inter_by_both1.txt.uu")

#s1 = readLines("../data/inter_by_both1.txt")
#s2 = readLines("../data/inter_by_both1.txt.old")

#s1 = readLines("../data/spca7_5pcs_inter2/PC5.txt")
#s2 = readLines("../data/spca7_15pcs_inter2/PC5.txt")

#s1 = readLines("../data/inter_by_both1.txt")
s1 = readLines("../data/marker_genes.txt")
s2 = readLines("../data/spca7_10pcs_inter1_sux_0.01_uu_0.6_0.7/PC10.txt")

#s1 = readLines("../data/spca7_10pcs_inter1_sux_0.01_cc_0.6_0.7/all_pcs.txt")
#s2 = readLines("../data/spca7_10pcs_inter1_sux_0.01_uu_0.6_0.7/all_pcs.txt")

inter = intersect(s1, s2)

print("length(s1)")
print(length(s1))
print("length(s2)")
print(length(s2))
print("length(inter)")
print(length(inter))

#writeLines(inter, "../data/spca7_10pcs_inter1_sux_0.01_cc_0.6_0.7/all_pcs_cc_uu.txt")
#writeLines(inter, "../data/marker_genes2.txt")
#writeLines(inter, "../data/spca7/inter_s1_s2.txt")
