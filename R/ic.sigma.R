`ic.sigma` <-
function(phy, matrix)
{
	f<-function(x) pic(x, phy)
	ic<-apply(matrix, 2, f)
	r<-crossprod(ic, ic)/nrow(ic)
	return(r)
}

