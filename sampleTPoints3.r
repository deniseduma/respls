sampleTPoints3 = function(apoints) {
	
	nOmit = 8
	
	j = 1
	allpoints = apoints
	omit = character(nOmit)
	while (j <= nOmit) {
		omit[j] = sample(allpoints, size=1)
		allpoints = setdiff(allpoints, omit) 
		j = j + 1
	}
	
	index=which(apoints %in% omit)

	cat("\n")
	print("omit")
	print(omit)
	print("index")
	print(index)
	
	return(index)
}	


