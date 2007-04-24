\name{delta.tree}
\alias{delta.tree}
\alias{lambda.tree}
\alias{tworate.tree}
\alias{linearchange.tree}
\alias{rescale.tree}
\title{Tree transformations}
\description{
	Apply various transformation to the branches of a phylogenetic tree;
	These are meant to correspond with changing the model of phenotypic evolution 
	for discrete or continuous characters
	}
\usage{
delta.tree(phy, delta, rescale = F)
lambda.tree(phy, lambda)
tworate.tree(phy, breakPoint, f) 
linearchange.tree(phy, endf)
rescale.tree(phy, total.length)
}
\arguments{
  \item{phy}{Phylogenetic tree}
  \item{delta}{Delta value.  This is the transformation suggested by Pagel as a test for a slow-down or speed-up.  Value of delta >1 disproporionately increase the length of external nodes (speed-up), where as values <1 disproportionately increase the length of internal nodes (slow-down).}
  \item{lambda}{Lambda value. This performs a transformation that was advocated by Pagel as a test of phylogenetic signal.  Values of lambda are alwasys less than one.  Each internal edge is multiplied by lambda, having the effect of reducing or eliminating phylogenetic signal.}
  \item{rescale}{if TRUE, rescale tree to a total length of one before applying transformation}
  \item{breakPoint}{This splits the tree into two parts determined by the breakpoint}
  \item{f}{Branch lengths to the right of the breakpoint (i.e., toward the tips) are multiplied by some factor f.  Thus, if f is greater than 1, these branches are made longer, whereas if f is less than 1 these branches are made shorter.}
  \item{endf}{This causes the rate of change to increase or decrease linearly, with the ratio of the beginning rate to the ending rate equal to endf}
  \item{total.length}{This is a general function that rescales the root to tip length of the tree to some pre-specified value.} 
}
\details{
}
\value{
	Phylogenetic tree in ape format
}
\references{PAGEL}
\author{Luke J. Harmon}
\note{ 

}
\seealso{}
\examples{


data(geospiza)
attach(geospiza)

delta.tree(geospiza.tree, 0.5)->g2
plot(g2)
lambda.tree(geospiza.tree, 0.5)->g3
plot(g3)
tworate.tree(geospiza.tree, 0.5, 0.5)->g4
plot(g4)
linearchange.tree(geospiza.tree, 0.1)->g5
plot(g5)

}
\keyword{ manip }
