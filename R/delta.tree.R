############################################################
###############TREE TRANSFORMATIONS#########################
############################################################

###LAMBDA###
lambda.tree <- function(phy, lambda)
{
	original.rtt <- max(branching.times(phy))
	ltree <- phy
	ltree.old <- new2old.phylo(phy) #needed due to change in ape tree format
	ltree$edge.length <- ltree$edge.length * lambda #shortens internal branches in proportion to lambda
	t <- original.rtt * (1-lambda) #needed to rescale tree back to its original root-to-tip length
	which(ltree.old$edge[,2] > 0) -> terminal.edges #identifies termal edges to be extended for recovery of original rtt length; based on old tree format in which nodes are coded as negative values
	ltree$edge.length[terminal.edges] + t -> ltree$edge.length[terminal.edges] #extends terminal edges
	return(ltree)
}

###DELTA###
delta.tree<-function(phy, delta, rescale=F)
{
	tmp<-as.numeric(phy$edge)
	times<-branching.times(phy)	
	max(times)->original.rtt.length
	times=max(times)-times
	max(times)
	res<-phy
	for(i in 1:length(phy$edge.length)) {
		bl<-phy$edge.length[i]
		age=times[which(names(times)==phy$edge[i,1])]
		res$edge.length[i]<-(age+bl)^delta-age^delta
	}
	if(rescale==T) res<-rescale.tree(res,original.rtt.length)
	res
}

###TWORATE#####
tworate.tree<-function(phy, breakPoint, f) 
{
	times<-branching.times(phy)	
	for(i in 1:length(phy$edge.length)) {
		bl<-phy$edge.length[i]
		age=times[which(names(times)==phy$edge[i,1])] #gets tip to node length
		if((age-bl)<breakPoint) #identifies branches that are on the tip side of the break-point (i.e., young)
			phy$edge.length[i]<-(age-min(age, breakPoint))*1+(min(age, breakPoint)-(age-bl))*f #If an edge is entirely to the right of the breakpoint it is simply multipled by f.  However, if the edge extends across the breakpoint, we want to leave the part to the left of the BP unchanged (this is the first part of this equation, which multiplies this part by 1).  We then want to multiply the part of the edge that is to the right of the BP by f (this is the second part of this line).
		}
	phy
}

###LINEAR CHANGE###
linearchange.tree<-function(phy, endf)
{
	times<-branching.times(phy)
	names(times)<-(as.numeric(names(times)))
	for(i in 1:length(phy$edge.length)) {
		bl<-phy$edge.length[i]
		age=times[which(names(times)==phy$edge[i,1])]
		mid=age-bl/2 
		rate=1+(endf-1)*(1-mid/max(times)) 
		phy$edge.length[i]<-phy$edge.length[i]*rate
	}
	phy	
}

###RESCALING###
##
rescale.tree<-function(phy, total.length)
{
	d<-max(branching.times(phy))
	phy$edge.length<-(phy$edge.length/d)*total.length
	phy
}

exponentialchange.tree<-function (phy, r) 
{
	
	if(r==0) return(phy)
    times <- branching.times(phy)
    names(times) <- (as.numeric(names(times)))
    for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age = times[which(names(times) == phy$edge[i, 1])]
        t1 = max(times) - age
        t2 = t1+bl
        phy$edge.length[i] = (exp(r*t2)-exp(r*t1))/(r)
    }
    phy
}

speciational.tree<-function(phy)
{
	phy$edge.length[]=1
	phy	
}