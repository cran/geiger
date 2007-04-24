`dtt` <-
function(phy, data, disp="avg.sq")
{
	phy2<-phy
	phy<-new2old.phylo(phy)
	
	result<-numeric()
	

	node.depth<-branching.times(phy2);
	stem.depth<-numeric();
	stem.depth[1]<-node.depth[1];
	for(i in 2:phy2$Nnode) {
		anc<-which(as.numeric(phy$edge[,2])==-i)
		stem.depth[i]<-node.depth[names(node.depth)==phy2$edge[anc,1]]
	}
		
	ltt<-sort(node.depth, decreasing=TRUE)
	node.depth<-node.depth/max(ltt);
	stem.depth<-stem.depth/max(ltt);
	ltt<-ltt/max(ltt);
	if(length(dim(data))==2) {
		d<-tip.disparity(phy2, data, disp);
		result[1]<-d[1]
		for(i in 2:length(ltt)) {
			x<-d[stem.depth>=ltt[i-1]&node.depth<ltt[i-1]]
			if(length(x)==0) result[i]=0
			else result[i]<-mean(x);
		}
		result[length(ltt)+1]<-0;
		if(result[1]>0)
			result<-result/result[1];
			
	} else {
		if(length(dim(data))!=3)
			stop("Error in data");
		
		for(i in 1:dim(data)[3]) {
			pp<-data[,,i]
			d<-tip.disparity(phy2, pp, disp);
			y<-numeric()
	
			y[1]<-d[1]
			for(j in 2:length(ltt)) {
				x<-d[stem.depth>=ltt[j-1]&node.depth<ltt[j-1]]
				if(length(x)==0) y[j]=0
				else y[j]<-mean(x);
			}
			y[length(ltt)+1]<-0;
			if(y[1]>0)
			y<-y/y[1];
			
			result<-cbind(result, y)
		}
	}
	
	return(result);	
}


