#!/usr/bin/Rscript

rm(list=ls())
setwd(".")
print(getwd())

library(MASS)
library(GeneCycle)
library(kernlab)
library(Matrix)
library(gplots)
library(multicore)
library(elasticnet)

args=commandArgs(TRUE)

outdir = "../data/spca7/"

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

#FIXME
sstring = "bs3|bs4|bs1_2|bs5"
regen_uu = regen_uu[, grep(sstring, colnames(regen_uu), invert=TRUE)]

geneSubset = readLines("../data/inter_by_both2.txt")

regen_cc = regen_cc[which(rownames(regen_cc) %in% geneSubset), ]
regen_uu = regen_uu[which(rownames(regen_uu) %in% geneSubset), ]
print("dim(regen_cc)")
print(dim(regen_cc))
print("dim(regen_uu)")
print(dim(regen_uu))

regen_big = cbind(regen_cc, regen_uu)
regen_big = t(regen_big)
#regen_big = log10(1 + regen_big)
print("dim(regen_big)")
print(dim(regen_big))

regen_t =  regen_big[which(rownames(regen_big) %in% c(tpoints_cc)), ]
regen_c =  regen_big[which(rownames(regen_big) %in% c(cpoints_cc)), ]

##Normalization of predictors
regen_t = scale(regen_t, center=TRUE, scale = FALSE)
regen_c = scale(regen_c, center=TRUE, scale = FALSE)
var_t = apply(regen_t, 2, function(l){var(l)})
var_c = apply(regen_c, 2, function(l){var(l)})
res_t = sort(var_t, decreasing=TRUE, index.return=TRUE)
res_c = sort(var_c, decreasing=TRUE, index.return=TRUE)
#regen_t = regen_t[, which(1:ncol(regen_t) %in% res_t$ix[1:2000])]
#regen_c = regen_c[, which(1:ncol(regen_c) %in% res_t$ix[1:2000])]
regen_t = regen_t[, res_t$ix[1:2000]]
regen_c = regen_c[, res_t$ix[1:2000]]
#Scale vars
#norm_t <- sqrt(drop(rep(1, nrow(regen_t)) %*% (regen_t^2)))
#regen_t <- scale(regen_t, FALSE, norm_t)
#norm_c <- sqrt(drop(rep(1, nrow(regen_c)) %*% (regen_c^2)))
#regen_c <- scale(regen_c, FALSE, norm_c)

cat("\n")
print("dim(regen_t)")
print(dim(regen_t))
print(rownames(regen_t))
print("dim(regen_c)")
print(dim(regen_c))
print(rownames(regen_c))
cat("\n")

timepts = rownames(regen_t)
genes = colnames(regen_t)
bs1_rows = grep("bs2", rownames(regen_t), invert=TRUE)
bs2_rows = grep("bs2", rownames(regen_t), invert=FALSE)

dim = dim(regen_t)
n = dim[1]
p = dim[2]

regen = regen_t[bs1_rows, ]
avgp(regen, angular=TRUE)
res = dominant.freqs(regen, m=3)
print("dim(res)")
print(dim(res))
print("any constant gene")
print(sum(is.constant(regen)))

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
	plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))

	# TODO: why this scaling is necessary?
	plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
        plot(plot.data, t="h", lwd=2, main="", xlab="Frequency (Hz)", ylab="Strength", xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}


plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
  Xk.h <- rep(0,length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
  points(ts, harmonic.trajectory, type="l", col=color)
}

#par(mfrow=c(4,4))
#ss = sample(1:p, 4)
#for (i in 1:4) {
#	traj1 = regen_t[bs1_rows, ss[i]]
#	traj2 = regen_t[bs2_rows, ss[i]]
#	traj3 = 0.5*traj1 + 0.5*traj2
#	X.k.1 = fft(traj1)
#	X.k.2 = fft(traj2)
#	X.k.3 = fft(traj3)
#	print("length(X.k.1)")
#	print(length(X.k.1))
#	plot(traj1, type="l")
#	plot(traj2, type="l")
	#plot(traj3, type="l")
#	plot.frequency.spectrum(X.k.1)
#	plot.frequency.spectrum(X.k.2)
	#plot.frequency.spectrum(X.k.3)
#}	


