\name{fit.discrete}
\alias{fit.discrete}
\title{ Model fitting for discrete comparative data }
\description{
	Fits macroevolutionary models to phylogenetic trees
}
\usage{
fit.discrete(tip.data, phy, lambda=FALSE, delta=FALSE, start=1.0)
}
\arguments{
  \item{tip.data}{ Character states for tip data; names should match, and data should be a factor with levels }
  \item{phy}{ Phylogenetic tree in ape format }
  \item{lambda}{ Estimate lambda? }
  \item{delta}{ Estimate delta? }
  \item{start}{ Starting guess for q }
}
\details{
}
\value{
 	Returns maximum likelihood value for q and selected parameters, along with the likelihood score
	Sometimes, parameters are confounded!
}
\references{ }
\author{ Luke J. Harmon and R. E. Glor }
\note{ 
}
\seealso{ }
\examples{
data(geospiza)
attach(geospiza)
drop.tip(geospiza.tree, "olivacea")->g.tree

gb<-as.factor(geospiza.data[,1]>4.2)
names(gb)<-rownames(geospiza.data)

fit.discrete(gb, g.tree)
fit.discrete(gb, g.tree, lambda=TRUE)
fit.discrete(gb, g.tree, delta=TRUE)


}
\keyword{ arith }
