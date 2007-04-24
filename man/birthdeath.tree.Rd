\name{birthdeath.tree}
\alias{birthdeath.tree}
\title{Birth-death tree simulator}
\description{
	Simulates phylogenetic tree under a uniform birth-death process	
}
\usage{
birthdeath.tree(b, d, time.stop = 0, taxa.stop = 0, seed = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{Per-lineage birth (speciation) rate}
  \item{d}{Per-lineage death (extinction) rate}
  \item{time.stop}{Stopping time}
  \item{taxa.stop}{Maximum number of taxa}
  \item{seed}{Random number seed; if seed=0 (default) then random number 
              generator is seeded based on the clock}
}
\details{
}
\value{
	Phylogenetic tree in ape format
}
\references{Geiger}
\author{Luke J. Harmon and Jason Wier}
\note{ 
}
\seealso{
}
\examples{

# Pure-birth tree
p1<-birthdeath.tree(b=0.1, d=0, time.stop=20)
plot(p1)

# Birth-death tree with extinct taxa
p2<-birthdeath.tree(b=0.2, d=0.05, time.stop=20)
plot(p2)

# Previous tree with extinct taxa removed
prune.extinct.taxa(p2)->p3
if(!is.null(p3)) plot(p3) else cat("No survivors\n")
#Sometimes there are no extant taxa

}
\keyword{arith}
\keyword{graphs}
