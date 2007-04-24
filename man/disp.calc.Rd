\name{disp.calc}
\alias{disp.calc}
\title{Calculate disparity}
\description{
	Calculate disparity of a set of points	
}
\usage{
disp.calc(data, disp = "avg.sq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{n x m matrix of data points, where n = number of points and m = number of dimensions}
  \item{disp}{disparity measure: avg.sq, avg.manhattan, nb.states}
}
\details{

}
\value{
	disparity
}
\references{Foote}
\author{Luke J. Harmon}
\note{

}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)

disp.calc(geospiza.data)

}

\keyword{ arith }
