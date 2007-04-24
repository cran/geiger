\name{dtt}
\alias{dtt}
\title{Disparity-through-time}
\description{
	Creates a disparity-through-time plot for either a single data set or multiple 
	data sets (for example, from simulations)
}
\usage{
dtt(phy, data, disp = "avg.sq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{Phylogenetic tree in 'phylo' format}
  \item{data}{Data - either actual or simulated}
  \item{disp}{Disparity metric - see disp.calc}
}
\details{
}
\value{
	Average disparity for clades whose stem crosses each time interval in the tree
}
\references{Harmon et al. 2003}
\author{ Luke J. Harmon }
\note{ 
}
\seealso{ }
\examples{

data(geospiza)
attach(geospiza)

drop.tip(geospiza.tree, "olivacea")->g.tree

disp.data<-dtt(g.tree, geospiza.data)
ltt<-sort(branching.times(g.tree), decr=TRUE)
ltt<-c(0, (max(ltt)-ltt)/max(ltt));
plot(ltt, disp.data, type="l", lwd=2)

}

\keyword{ arith }
