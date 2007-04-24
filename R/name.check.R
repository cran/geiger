`name.check` <-
function(x, phy)
{
	t<-phy$tip.label;
	d<-rownames(x);
	r1<-t[is.na(match(t,d))]
	r2<-d[is.na(match(d,t))]
	r<-list(sort(r1), sort(r2))
	names(r)<-cbind("Tree.not.data", "Data.not.tree")
	if(length(r1)==0 && length(r2)==0) return("OK")
	else return(r)
}

