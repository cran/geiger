`tip.disparity` <-
function(phy, data, disp="avg.sq")
{
	if (class(phy) != "phylo")
		stop("object \"phy\" is not of class \"phylo\"");
		
	nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    result<-numeric();
    
    for(i in 1:nb.node) {
    	l<-node.leaves(phy, nb.tip+i);
    	d<-data[match(l, row.names(data)),];
    	result[i]<-disp.calc(d, disp);
    }
    return(result);	  	
}

