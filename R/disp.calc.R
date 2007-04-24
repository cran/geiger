`disp.calc` <-
function(data, disp="avg.sq")
{
	if(disp=="avg.sq") {
		d<-dist(data, method="euclidean")^2
		r<-mean(d)
	}
	else if(disp=="avg.manhattan") {
		d<-dist(data, method="manhattan")
		r<-mean(d)
	}
	else if(disp=="num.states") {
		f<-function(x) length(unique(x))
		d<-apply(data, 2, f)
		r<-mean(d)
	}
	else r<-0;
	return(r)
}

