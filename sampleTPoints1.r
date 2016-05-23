
sampleTPoints = function(apoints, nOmit) {

	counts = list()
	counts[c("0", "24", "48", "54", "60", "66", "72", "96", "120", "144", "168")] = c(5, 6, 4, 2, 2, 2, 5, 5, 6, 4, 5)
	myl = list()
	myl[c("0", "24", "48", "54", "60", "66", "72", "96", "120", "144", "168")] = rep(0, 11)
	
	j = 1
	allpoints = apoints
	omit = character(nOmit)
	while (j <= nOmit) {
		mys = sample(allpoints, size=1)
		
		mys2 = mys
		if (grepl("cc_", mys) == TRUE) {
			mysplt = strsplit(mys, "cc_")
			mys2 = mysplt[[1]][2]
		}
		mysplt = strsplit(mys2, type, fixed = TRUE)
		ss = mysplt[[1]][1]  
		
		while ((myl[[ss]] + 1) >= counts[[ss]]) {
			mys = sample(allpoints, size=1)
			
			mys2 = mys
			if (grepl("cc_", mys) == TRUE) {
				mysplt = strsplit(mys, "cc_")
				mys2 = mysplt[[1]][2]
			}
			mysplt = strsplit(mys2, type, fixed = TRUE)
			ss = mysplt[[1]][1]  
		}
		
		j = j + 1
		omit[j] = mys
		myl[[ss]] = myl[[ss]] + 1
		allpoints = setdiff(allpoints, omit) 
	}
	
	print("omit")
	print(t(omit))
	print("counts")
	print(t(counts))
	print("myl")
	print(t(myl))
	
	omit

}	


