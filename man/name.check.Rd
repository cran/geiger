\name{name.check}
\alias{name.check}
\title{ Compares taxa in data and tree }
\description{
Checks for taxa present in data set but not in tree, and vice-versa
}
\usage{
name.check(x, phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Data }
  \item{phy}{ Tree }
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

name.check(geospiza.data, geospiza.tree)

}

\keyword{ arith }
