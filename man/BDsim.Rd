\name{BDsim}
\alias{BDsim}
\title{Birth-death population simulator}
\description{
	Simulates species richness (or population growth) under a uniform birth-death process	
}
\usage{
BDsim(nStart, b, d, times) 
}
\arguments{
  \item{nStart}{Number of taxa at starting time zero}
  \item{b}{Per-lineage birth (speciation) rate}
  \item{d}{Per-lineage death (extinction) rate}
  \item{times}{Vector of times where extant species are counted}
}
\details{
}
\value{
	Population size at each time in TIMES
	}
\references{Geiger}
\author{Richard E. Glor and Luke J. Harmon}
\note{ 
}
\seealso{
}
\examples{

pop1=BDsim(nStart=10, b=0.1, d=0, times=1:10)
pop2=BDsim(nStart=10, b=0, d=0.1, times=1:10)
pop3=BDsim(nStart=10, b=0.1, d=0.1, times=1:10)

plot(pop1, type="l", ylim=c(0,20))
lines(pop2, col="red")
lines(pop3, col="blue")

}
\keyword{arith}
