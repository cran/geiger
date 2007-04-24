\name{area.between.curves}
\alias{area.between.curves}
\title{Area between two curves}
\description{
	Finds the area between two curves up to a cutoff point
	(used to calculate the MDI statistic from Harmon et al. 2003
}
\usage{
area.between.curves(x, f1, f2, xrange = c(0,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of values for the x axis}
  \item{f1}{Y values for curve 1}
  \item{f2}{Y values for curve 2}
  \item{xrange}{Range of x values over which area is calculated}
}
\details{
	Area is calculated between the two curves f1 and f2 for x values between 0 and cutoff
	If f1 is greater than f2, areas are positive; otherwise, areas are negative
}
\value{
	Area between f1 and f2 between 0 and cutoff
}
\references{Harmon et al. 2003}
\author{Luke J. Harmon}
\note{
}
\examples{

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

area.between.curves(ltt, mean.sims, disp.data, xrange=c(0, 2/3))
}

\keyword{arith}
