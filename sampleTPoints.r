sampleTPoints = function(apoints) {
	
	count = 0
	index = numeric(length(apoints)/2)
	if (length(grep('cc', apoints[1]))) {
		tpoints = 16
	} else {
		tpoints = 22
	}
	for (j in seq(1,tpoints,2)) {
		count = count + 1
		index[count] = sample(c(j, j+1), size=1)
	}
	return(index)
}	

