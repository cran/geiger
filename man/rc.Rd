\name{rc}
\alias{rc}
\title{ Relative cladogenesis test }
\description{
	Carries out relative cladogenesis test for all slices through the tree.	
}
\usage{
rc(phy, plot = T)
}
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{plot}{ Plot tree with significant branches highlighted? }
}
\details{
}
\value{
  Table of results with three columns: number of ancestors, maximum descendents, p-value
}
\references{ XXX }
\author{ Luke J. Harmon }
\note{ 
}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)

rc(geospiza.tree)

}

\keyword{ arith }
