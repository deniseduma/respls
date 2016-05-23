#!/usr/bin/Rscript


library(igraph)

# your code
#require(igraph)
#g1 <- erdos.renyi.game(10, 0.5)
#V(g1)$time <- 1:10
# adding a weight (difference in times)
#E(g1)$weight <- apply(get.edgelist(g1), 1, function(i){abs(V(g1)$time[i[1]]-V(g1)$time[i[2]])})
# calculate adjacency and/or shortest path. 
#get.adjacency(g1, attr="weight")
#shortest.paths(g1)

new = 0
p = 500
outdir = "test/" #"sim3/"

## Create scale-free graph
#algorithm = c("psumtree", "psumtree-multiple", "bag")
#adj_mat = matrix(c(0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0),5,5)
#g_adj_mat = graph_from_adjacency_matrix(adj_mat, mode = "undirected")
if (new) {
	g = sample_pa(p, power = 1.3, m = 1, out.dist = NULL, out.seq = NULL,out.pref = FALSE, zero.appeal = 1, directed = FALSE, algorithm = "psumtree", start.graph = NULL)
	## Save graph to file
	write_graph(g, paste(outdir, "igraph.txt",sep=""), format = "edgelist")
} else {
	## Read graph from file
	g = read_graph(paste(outdir,"igraph.txt",sep=""), format="edgelist", directed=FALSE)
}
#print("g")
#print(str(g))
#degree_distribution(g)
print("degree_distribution(g)")
ff=degree_distribution(g)
plot(log(1:length(ff)), log(ff))
print(paste("R^2 ", cor(log(1:length(ff)), log(ff), use="pairwise.complete.obs"), sep=""))

## Vizualize it
if (new) {
	pdf(paste(outdir, "igraph.pdf",sep=""))
	plot(g,vertex.size=8,vertex.label=NA,layout=layout.fruchterman.reingold)
	dev.off()
}	

## Communities
if (new) {
	#cls = cluster_label_prop(g)
	#cls = cluster_walktrap(g)
	#cls = cluster_leading_eigen(g, steps=4)
	#cls = cluster_fast_greedy(g)
	cls = cluster_edge_betweenness(g) 
	# Save communities to file
	saveRDS(cls, paste(outdir,"communities.txt",sep=""))
} else {
	cls = readRDS(paste(outdir, "communities.txt", sep=""))
}

ncom = length(cls)
print(paste("length(cls) ", ncom, sep=""))
sz = sizes(cls)
print("sizes(cls)")
print(sz)
cls2 = membership(cls)
#cmts = communities(cls)

## Smaller no of communities 
ncom=10
cls2 = cut_at(cls, no=ncom)
sz = sapply(1:ncom, function(j){sum(cls2==j)})
print("sizes(cls2)")
print(sz)

if (new) {
	sz = as.vector(sz)
	write.table(sz, file = paste(outdir,"c_sizes.txt",sep=""), quote = F,row.names=F,col.names=F)
}

if (new) {
	## Plot communities
	pdf(paste(outdir,"communities.pdf",sep=""))
	##plot(cls,g)
	colors = rainbow(ncom)
	plot(g,vertex.color=colors[cls2],vertex.size=8,vertex.label=NA,layout=layout.fruchterman.reingold)
	dev.off()
}

if (new) {
	## Compute V
	V = matrix(0, p, ncom)
	for (j in 1:ncom) {
		idx = which(cls2==j)
		V[ idx, j] = 1	
	}
	# Write V to file
	write.table(V, file = paste(outdir,"V.txt",sep=""),quote = F,row.names=F,col.names=F)
} else {
	# Read V from file
	V = read.table(paste(outdir,"V.txt",sep=""),header=FALSE,quote = NULL)
}	
print("V")
print(dim(V))

if (new) {
	## Adj matrix & Laplacian
	A = as.matrix(get.adjacency(g))
	L = laplacian_matrix(g, normalized = FALSE, weights = NULL,sparse = FALSE)
	# Write A & L
	write.table(A, file = paste(outdir, "A.txt",sep=""),quote = F,row.names=F,col.names=F)
	write.table(L, file = paste(outdir,"L.txt",sep=""),quote = F,row.names = F,col.names = F)
} else {
	# Read A & L
	A = as.matrix(read.table(paste(outdir, "A.txt",sep=""), header=FALSE, quote = NULL))
	L = as.matrix(read.table(paste(outdir, "L.txt",sep=""), header=FALSE, quote = NULL))
}
print("dim(A)")
print(dim(A))
print("dim(L)")
print(dim(L))

try(chol(L))


