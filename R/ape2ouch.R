`ape2ouch` <-
function(phy, data)
{
	phy<-new2old.phylo(phy)
	n<-length(phy$edge.length)+1
	node<-1:n
	species<-character(n)
	species[1]<-NA
	species[which(as.numeric(phy$edge[,2])<0)+1]<-NA
	species[which(as.numeric(phy$edge[,2])>0)+1]<-phy$tip.label
	ancestor<-match(phy$edge[,1], phy$edge[,2])
	ancestor[is.na(ancestor)]<-0
	ancestor<-ancestor+1
	ancestor<-c(NA, ancestor)
	time<-numeric(n)
	time[1]<-0.0;
	for(i in 2:length(time))
		time[i]<-time[ancestor[i]]+phy$edge.length[i-1];
	d<-rep(NA, length(species))
    d[match(names(data),species)]<-data

	obj<-list(d=d, node=node, species=species, ancestor=ancestor, time=time)
	obj
}



