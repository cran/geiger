##evenQ is an internal function of GEIGER
##This function makes a template for the calculate of a rate matrix that will set all transitions to the same value (based on the total number of states).  One can then multiply the resulting matrix by the overall rate to get the rate matrix for a particular analysis.


`evenQ` <-
function(n)

{

	q<--diag(n)

	q[q==0]<-1/(n-1)

	return(q)

}

