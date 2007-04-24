\name{geiger-package}
\alias{geiger-package}
\alias{geiger}
\docType{package}
\title{
GEIGER
}
\description{
 A package for macroevolutionary simulation and estimating parameters 
 related to diversification from comparative phylogenetic data.
}
\details{
\tabular{ll}{
Package: \tab geiger\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-01-19\cr
License: \tab GPL version 2 or greater?\cr
}

}
\author{
Luke J. Harmon, Jason Wier, Chad Brock, and Wendell Challenger

Maintainer: Luke Harmon <harmon@zoology.ubc.ca>
}
\references{
}
\seealso{

}
\examples{

# Disparity-through-time, as in Harmon et al. 2003
data(geospiza)
attach(geospiza)

drop.tip(geospiza.tree, "olivacea")->g.tree
disp.data<-dtt(g.tree, geospiza.data)
s<-ic.sigma(g.tree, geospiza.data)
sims<-sim.char(g.tree, s, 100)

disp.sims<-dtt(g.tree, sims)
mean.sims<-rowSums(disp.sims)/ncol(disp.sims)

ltt<-sort(branching.times(g.tree), decr=TRUE)
ltt<-c(0, (max(ltt)-ltt)/max(ltt));
plot(ltt, disp.data)
lines(ltt, disp.data)
lines(ltt, mean.sims)


# Simulation tester - verifies that simulations produce expected vcv structure
sims<-sim.char(geospiza.tree, as.matrix(1.0), n=100000)
m<-var(t(sims[,1,]))
m2<-vcv.phylo(geospiza.tree)
round(m-m2)
# Equal to zero because the vcv matrices are as expected.

#Interface with ouch package
data(geospiza)
attach(geospiza)
name.check(geospiza.data, geospiza.tree)->r
drop.tip(geospiza.tree, r[[1]])->g.tree

data<-geospiza.data[,1]
names(data)<-rownames(geospiza.data)

ape2ouch(g.tree, data)->geospiza.ouch

#Interface with ouch package
library(ouch)

brown.fit(geospiza.ouch$d, geospiza.ouch$node, geospiza.ouch$ancestor, geospiza.ouch$time)

# Reconstruct evolutionary vcv matrix
v<-cbind(c(4.0, -1.5),c(-1.5, 9.0))
n<-1000
sims<-sim.char(geospiza.tree, v, n)
r<-array(dim=c(2, 2, n))
for(i in 1:n)
	r[,,i]<-ic.sigma(geospiza.tree, sims[,,i])
v.recon<-apply(r, c(1,2), mean)
v-v.recon

# Discrete character evolution
q<-list(rbind(c(-.01, .01), c(.01, -.01)))
sims<-sim.char(geospiza.tree, q, model="discrete", n=100)

}
\keyword{arith}
\keyword{graphs}