###Felsenstein's pruning algorithm
likelihood.discrete<-function(phy, tip.data, q, delta=1, lambda=1,  kappa=1, endRate=1, linear=F, breakPoint=0, f=1, rtt.rescale=0, total.rescale=F, returnFull=F)
{
	
	if(!is.factor(tip.data)) tip.data<-factor(tip.data)
	Q<-evenQ(nlevels(tip.data))*q
	if (class(phy) != "phylo")
		stop("object \"phy\" is not of class \"phylo\"");
	#new2old.phylo(phy)->phy ##converts back to old ape tree format with negative values denoting internal nodes
	tmp <- as.numeric(phy$edge)
	nb.tip <- max(tmp) #number of tips
	nb.node <- -min(tmp) #number of internal nodes
	nb.states <- nlevels(tip.data) #numbers of states in character
	l <- matrix(0, nrow=nrow(phy$edge), ncol=nb.states) #makes matrix that will store likelihoods of each character state at each node
	root <- numeric(nb.states) #makes vector that will store root likelihoods
	m <- match(phy$tip.label, names(tip.data)) ##identifies elements of tip.data matrix that corresponds with the tip.label
	if (delta != 1)
		deltaTree(phy, delta) -> phy;
	if (lambda != 1)
		lambdaTree(phy, lambda) -> phy;
	if (kappa != 1)
		kappaTree(phy, kappa) -> phy;
	if (endRate != 1) {
		if(breakPoint!=0) {
			tworateTree(phy, breakPoint, endRate) -> phy;
		} else if(linear==T) {
			linearchangeTree(phy, endRate) -> phy;
		} else exponentialchangeTree(phy, endRate)->phy;
	}


	
	#When comparing deltas across different qs, it might be useful to rescale the total tree length to one	
	if(rtt.rescale!=0)	
		rescaleTree(phy, rtt.rescale) -> phy
		
	new2old.phylo(phy)->phy
	for(i in 1:nrow(phy$edge)) #for each edge
		if(as.numeric(phy$edge[i,2])>0) l[i,tip.data[m[as.numeric(phy$edge[i,2])]]] <- 1.0 #if the edge is connected to a terminal taxon, you set the likelihood of the tip value equal to 1 and all others equal to zero.
		times <- branching.times(old2new.phylo(phy)) #get node to tip distances
		-1*(1:(max(as.numeric(names(times)))-min(as.numeric(names(times)))+1))->names(times)
		times = max(times) - times #convert into root to node tips
		if(total.rescale) {
			sum(phy$edge.length) -> total.tree
			phy$edge.length <- phy$edge.length/total.tree
		}
	while(1) {
		
		if(sum(as.numeric(phy$edge[,2])>0)==2) break 
	
		#obtain ancestors of current tips
		x <- match(phy$edge[,1], phy$edge[,2])[as.numeric(phy$edge[,2])>0] #finds nodes connected to terminal taxa
		#find last node with two tip descendent
		a <- max(x[duplicated(x)])
		t <- which(phy$edge[,1]==phy$edge[a,2])
		bl <- phy$edge.length[t]
		age = times[which(names(times)==phy$edge[a,2])]
		l[a,] <- frag.like(l[t,], bl, Q)

		#next line effectively prunes out the tips just used
		phy$edge[a,2]<-1
		phy$edge[t,2]<-0

	
	}
	t <- which(as.numeric(phy$edge[,2])>0)
	bl <- phy$edge.length[t]
	root <- frag.like(l[t,], bl, Q)
	neglnl=-log(sum(root/nb.states))
	if(returnFull==F) {
		return(neglnl)
	} else return(list(neglnl=neglnl, root=root, l=l))
}
