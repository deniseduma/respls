#!/usr/bin/Rscript

rm(list=ls()) 
setwd("../data/nscumcomp_t/genes/")
print(getwd())

args=commandArgs(TRUE)
genes = readLines(args[1])
print("length(genes)")
print(length(genes))
print("length(unique(genes))")
print(length(unique(genes)))

