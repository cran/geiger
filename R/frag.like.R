#This function calculates the likelihood on one branch of a tree

`frag.like` <-
function(tip.like, bl, q)

{

	nb.states<-ncol(tip.like)

	r<-rep(1, nb.states)

	d<-length(bl)

	p<-list(d)

	for(i in 1:d)

		p[[i]]<-MatrixExp.simple(q*bl[i])

	for(i in 1:nb.states)

		for(j in 1:d) 

			r[i]<-r[i]*sum(p[[j]][i,]*tip.like[j,])

	return(r)

}

