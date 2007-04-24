\name{ape2ouch}
\alias{ape2ouch}

\title{Convert between ape and ouch formats}
\description{
  Converts phylogenetic trees between ape format (used by geiger, ape, and laser packages) and 
  ouch format (used by the ouch package)
}
\usage{
ape2ouch(phy, data)
}
\arguments{
  \item{phy}{Phylogenetic tree in ape format}
  \item{data}{Data vector with names that match tip species}
  
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
\references{Geiger2}
\author{Luke J. Harmon}

\examples{
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

}

\keyword{manip}



