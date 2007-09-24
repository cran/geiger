\name{name.check}
\alias{name.check}
\title{ Compares taxa in data and tree }
\description{
This function is a general tool for checking for concordance between a data file and a phylogenetic
tree.  For the data, names can be specified as the names of objects in the vector, rownames of the data 
array or as 'x.names'. The function finds and lists all taxa present in data set but not in tree, and vice-versa.
}
\usage{
name.check(phy, data, data.names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ an object of class "phylo" }
  \item{data}{ data for tips of the tree }
  \item{data.names} { names of the tips in the order of the data; if this is not given, names 
  					will be taken from the names or rownames of the object data}
}
\details{
}
\value{
  \item{Tree.not.data}{Taxa in tree but not data}
  \item{Data.not.tree}{Taxa in data but not tree}
  ...
}
\references{  }
\author{ Luke J. Harmon }
\note{ 

}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)

name.check(geospiza.tree, geospiza.data)

}

\keyword{ arith }
