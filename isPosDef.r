isPosDef<- function(M) { if ( all(M == t(M) ) ) {  # first test symmetric-ity
	if (  all(eigen(M)$values>  0) ) {TRUE}
		else {FALSE} 
	} else {FALSE}  # not symmetric
}
