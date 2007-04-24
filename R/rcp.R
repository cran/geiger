`rcp` <-
function(ni, n, k) 
{
	max<-floor((n-k)/(ni-1))
	sum=0
	for(v in 0:max) {
		term=(-1)^v*choose(k, v)*choose(n-v*(ni-1)-1, k-1)
		sum=sum+term
 	}
	return(1-(sum/choose(n-1, k-1)))
}

