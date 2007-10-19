\name{ape2ouch}
\alias{ape2ouch}

\title{Convert between ape and ouch formats}
\description{
  Converts phylogenetic trees between ape 'phylo' format (used by geiger, ape, and laser packages) and 
  ouch format (used by the ouch package).  The OUCH package can be used for testing the Ornstein-Uhlenbeck
  model of evolution, among other things.  If the tree and data file have taxa that don't match, they will
  be removed from the analysis.
}
\usage{
ape2ouch(phy, data, data.names=NULL)
}
\arguments{
  \item{phy}{ An object of class phylo}
  \item{data}{ Data vector }
  \item{data.names}{ Tip names for data vector that match tree species; ignored if data includes names}
  
}

\value{
  \code{ape2ouch} returns the phylogentic tree in ouch format as a list with the
  following components:
  \item{d}{Phenotypic values in ouch format}
  \item{node}{Node list for a phylogenetic tree, ouch format}
  \item{ancestor}{Ancestor list for a phylogenetic tree, ouch format}
  \item{time}{Branch times for a phylogenetic tree, ouch format}
  \item{species}{Species names for a phylogenetic tree, ouch format}
}
\references{ Butler, M.A. and A.A. King, 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. American Naturalist 164:683-695.}
\author{Luke J. Harmon}

\examples{
data(geospiza)
attach(geospiza)

ape2ouch(geospiza.tree, geospiza.data[,1], rownames(geospiza.data))->geospiza.ouch

#Interface with ouch package
library(ouch)

brown.fit(geospiza.ouch$d, geospiza.ouch$node, geospiza.ouch$ancestor, geospiza.ouch$time)

}

\keyword{manip}



