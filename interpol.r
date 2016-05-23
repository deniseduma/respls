#!/usr/bin/Rscript

rm(list=ls())
setwd(".")
print(getwd())

source("ss6.r")
source("spca6.r")
source("sampleTPoints3.r")

library(MASS)
library(glmnet)
library(Matrix)
library(gplots)
library(multicore)
library(elasticnet)
#library(ConsensusClusterPlus)

args=commandArgs(TRUE)
#type = "c"#args[1]
myncomp = 6#as.integer(args[2])
nOmit = 6 #as.integer(args[3])
noRuns = 1 #as.integer(args[4])
noSteps = 25
threshold = 0.8#threshold for stability selection 

#Read in input data and do some pre-processing
dir_t = paste("nscumcomp_", "t", sep="")
dir_c = paste("nscumcomp_", "c", sep="")

##Measure execution time
start = proc.time()

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

geneSubset = readLines("../data/inter_by_both.txt")
#regen_cc = regen_cc[geneSubset, ]
#regen_uu = regen_uu[geneSubset, ]

regen_cc = regen_cc[geneSubset, ]
regen_uu = regen_uu[geneSubset, ]
print("dim(regen_cc)")
print(dim(regen_cc))
print("dim(regen_uu)")
print(dim(regen_uu))

regen_big = cbind(regen_cc, regen_uu)
regen_big = t(regen_big)
print("dim(regen_big)")
print(dim(regen_big))

genes = colnames(regen_big)

regen_t =  regen_big[which(rownames(regen_big) %in% c(tpoints_cc)), ]
regen_c =  regen_big[which(rownames(regen_big) %in% c(cpoints_cc)), ]

dim = dim(regen_t)
n = dim[1]
p = dim[2]

#regen = log2(1+regen)
##Normalization of predictors
regen_t = scale(regen_t, center=TRUE, scale = FALSE)
regen_c = scale(regen_c, center=TRUE, scale = FALSE)

print("dim(regen_t)")
print(dim(regen_t))
#print(rownames(regen_t))
print("dim(regen_c)")
print(dim(regen_c))
#print((rownames(regen_c)))

tp = grep("bs2", rownames(regen_t))
ntp = setdiff(1:n, tp)
print("tp")
print(rownames(regen_t)[tp])
print("ntp")
print(rownames(regen_t)[ntp])

new_regen_t = cbind(regen_t[ntp, ], regen_t[tp, ])
new_regen_c = cbind(regen_c[ntp, ], regen_c[tp, ])
print("dim(new_regen_t)")
print(dim(new_regen_t))
print("dim(new_regen_c)")
print(dim(new_regen_c))

inter_regen_t = apply(new_regen_t, 2, function(col) {res = approx(1:8, y = col, method = "linear", n = 100, rule = 1, f = 0, ties = mean); return(res$y)})
inter_regen_c = apply(new_regen_c, 2, function(col) {res=approx(1:8, y = col, method = "linear", n = 100, rule = 1, f = 0, ties = mean); return(res$y)})
print("dim(inter_regen_t)")
print(dim(inter_regen_t))
print("dim(inter_regen_c)")
print(dim(inter_regen_c))

par(mfrow=c(4,2))
ss = sample(1:p, 8)
for (i in 1:8) {
	matplot(1:100, cbind(inter_regen_t[, ss[i]], inter_regen_t[, p + ss[[i]]]), type="l", xaxt="n", pch=16, xlab = "Time point", ylab = "Gene expression")
	axis(1, at=1:100)
}

