phy.anova<-function(y, factor, phy, nsim=1000)
{
	s<-mean(pic(y, phy)^2)
	
	a<-anova(lm(y~factor))
	f.data<-a[1,4]
	
	sims<-sim.char(phy, as.matrix(s), nsim=nsim)
	
	foo<-function(xx) anova(lm(xx~factor))[1,4]
	
	f.null<-apply(sims, 3, foo)

	cat("Standard ANOVA:\n")	
	print(a)
	
	cat("\n\nPhylogenetic p-value: \t")
	cat(sum(f.data>f.null)/(nsim+1))

		
}

phy.manova<-function(y, factor, phy, nsim=1000, test="Wilks")
{
	s<-ic.sigma(phy, y)
	
	m<-summary.manova(manova(as.matrix(y)~f), test=test)
	
	w.data<-m[[4]][1,2]
	
	sims<-sim.char(phy, s, nsim=nsim)
	
	foo<-function(xx) summary.manova(manova(as.matrix(xx)~f), test=test)[[4]][1,2]
	
	w.null<-apply(sims, 3, foo)

	cat("Standard MANOVA:\n")	
	print(m)
	
	cat("\n\nPhylogenetic p-value: \t")
	cat(sum(w.data>w.null)/(nsim+1))
	
}