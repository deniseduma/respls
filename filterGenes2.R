rm(list=ls())
setwd("/home/duma/HairCell/haircell")

source("samr.morefuns.R")

library(samr)

type = "cc"

## read cc
regen = read.table(paste("../data/", type, "_regen.csv", sep = ""), sep=",", header=T, row.names=1)
apoints = readLines(paste("../data/", type, "_timepoints.txt", sep=""))
cpoints = readLines(paste("../data/", type, "_timepoints_c.txt", sep=""))
tpoints = readLines(paste("../data/", type, "_timepoints_t.txt", sep=""))
colnames(regen) = apoints

print("dim(regen)")
print(dim(regen))

if (identical(type, "cc")) {
	sstring = "bs2"
} else if (identical(type, "uu")) {	
	sstring = "bs2|bs3|bs4|bs1_2|bs5"
}	

control = regen[, cpoints]
control_bs1 = control[, grep(sstring, cpoints, invert=TRUE)]
control_bs2 = control[, grep("bs2", cpoints, invert=FALSE)]
print("dim(control)")
print(dim(control))

treated = regen[, tpoints]
treated_bs1 = treated[, grep(sstring, tpoints, invert=TRUE)]
treated_bs2 = treated[, grep("bs2", tpoints, invert=FALSE)]
print("dim(treated)")
print(dim(treated))

regen = as.matrix(cbind(control_bs1, control_bs2, treated_bs1, treated_bs2))
print("dim(regen)")
print(dim(regen))

timepts = colnames(regen)
genes = rownames(regen)
ntimepts = length(timepts)/4
print("timepts")
print(timepts)

if (identical(type, "cc")) {
	#y=paste(c(rep(-1, 8), rep(1, 8)), "Time", c(1:8, 1:8), sep="")
	#start=c(1,9)
	#for(i in start){
	#	y[i]=paste(y[i],"Start",sep="")
	#}
	#end=c(8, 16)
	#for(i in end){
	#	y[i]=paste(y[i],"End",sep="")
	#}

	y=paste(c(rep(1, ntimepts), rep(1, ntimepts), rep(2, ntimepts), rep(2, ntimepts)), "Time", rep(1:ntimepts, 4), sep="")
	start=c(1, ntimepts + 1, 2*ntimepts + 1, 3*ntimepts + 1)
	for(i in start){
		y[i]=paste(y[i],"Start",sep="")
	}
	end=c(ntimepts, 2*ntimepts, 3*ntimepts, 4*ntimepts)
	for(i in end){
		y[i]=paste(y[i],"End",sep="")
	}
	
	#y=paste(c(rep(1,15),rep(2,15)),"Time",rep(c(1,2,3,4,5,1.1,2.5, 3.7, 4.1,5.5),3),sep="")
	#start=c(1,6,11,16,21,26)
	#for(i in start){
	#y[i]=paste(y[i],"Start",sep="")
	#}
	#for(i in start+4){
	#y[i]=paste(y[i],"End",sep="")
	#}

} else if (identical(type, "uu")) {
	#y = as.factor(c(seq(1, 30), seq(-1, -30)), "Time", rep(c(rep(1, 3), rep(2,4), rep(3, 2), rep(4,2), rep(5,2), rep(6,2), rep(7,3), rep(8, 3), rep(9, 4), rep(10, 2), rep(11,3)), 2))
	y=paste(c(rep(1, ntimepts), rep(1, ntimepts), rep(2, ntimepts), rep(2, ntimepts)), "Time", rep(1:ntimepts, 4), sep="")
	start=c(1, ntimepts + 1, 2*ntimepts + 1, 3*ntimepts + 1)
	for(i in start){
		y[i]=paste(y[i],"Start",sep="")
	}
	end=c(ntimepts, 2*ntimepts, 3*ntimepts, 4*ntimepts)
	for(i in end){
		y[i]=paste(y[i],"End",sep="")
	}

}

print("length(y)")
print(length(y))
print("y")
print(y)

#x = matrix(rnorm(1000*30),ncol=30)
#parse.time.labels.and.summarize.data(regen, y, "Two class paired timecourse", "slope")

samfit = SAM(regen, y, resp.type="Two class unpaired timecourse", genenames = genes,  nperms=1000, time.summary.type="slope", fdr.output=.05)
print("samfit$siggenes.table$genes.lo")
print(dim(samfit$siggenes.table$genes.lo))
print("samfit$siggenes.table$genes.up")
print(dim(samfit$siggenes.table$genes.up))
plot(samfit)

genes = c(samfit$siggenes.table$genes.lo[, 1], samfit$siggenes.table$genes.up[, 1])
print("length(genes)")
print(length(genes))

writeLines(genes, paste("../data/", type, "_by_both2.txt", sep=""))

