`rc` <-
function(phy, plot=T)
{

	
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
    
	nd<-branching.times(phy);
	node.depth<-c(nd, rep(0, nb.tip))
	names(node.depth)<-c(names(nd), as.character(1:nb.tip))
	p<-numeric(nb.node)
	max.desc<-numeric(nb.node)
	num.anc<-numeric(nb.node)
	node.name<-character(nb.node)

	stem.depth<-numeric();
	stem.depth[1]<-node.depth[1];
	
	for(i in 2:length(node.depth)) {
			anc<-which(phy$edge[,2]==names(node.depth)[i])
			stem.depth[i]<-node.depth[names(node.depth)==phy$edge[anc,1]]
	}
	

		
	ltt<-sort(nd, decreasing=TRUE)
	
	for(i in 2:length(ltt)) {
		nn<-stem.depth>=ltt[i-1]&node.depth<ltt[i-1]
		anc<-sum(nn)
		desc<-numeric(anc)
		pp<-numeric(anc)
		num.anc[i]<-anc
		for(j in 1:anc) {
			desc[j]<-length(node.leaves(phy, names(nn)[nn][j]))
		}
		max.desc[i]<-max(desc)
		p[i]<-rcp(max.desc[i], nb.tip, anc)
		node.name[i]<-names(ltt[i])
	}
	num.anc[1]<-1
	max.desc[1]<-nb.tip
	p[1]<-1
	res<-cbind(num.anc, max.desc, p)
	rownames(res)<-c("root", node.name[2:nb.node])
	
	if(plot) {
		labels<-character(length(phy$edge.length))
		names(labels)<-as.character(1:length(labels)+nb.tip)
		s<-which(res[,3]<0.05)
		labels[names(s)]<-"*"
		phy$node.label<-labels
		plot.phylo(phy, show.node.label=T)
	}
	return(res)
}

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


