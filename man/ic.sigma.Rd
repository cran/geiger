\name{ic.sigma}
\alias{ic.sigma}
\title{Estimate Evolutionary VCV Matrix}
\description{
	Uses independent contrasts to estimate the evolutionary variance-covariance matrix for 
	n quantitative characters
}
\usage{
ic.sigma(phy, matrix)
}
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{matrix}{ Matrix of continuous character data }
}
\details{
}
\value{
	Variance-covariance matrix
}
\references{ Revell et al. 2006 }
\author{ Luke J. Harmon }
\note{ 

}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)

drop.tip(geospiza.tree, "olivacea")->g.tree
ic.sigma(g.tree, geospiza.data)

}

\keyword{ arith }
