\name{sim.char}
\alias{sim.char}
\title{ Simulate character evolution }
\description{
  Simulates evolution of discrete or continuous characters on a phylogenetic tree
}
\usage{
sim.char(phy, model.matrix, nsims = 1, model = "brownian", root.state = 1)
}
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{model.matrix}{ Matrix describing model: either vcv matrix or q matrix }
  \item{nsims}{ Number of simulations to run }
  \item{model}{ "brownian", "speciational", or "discrete" }
  \item{root.state}{ Starting state at root }
}
\details{
}
\value{
 	Array of simulated data
}
\references{ }
\author{ Luke J, Harmon }
\note{ 
}
\seealso{  }
\examples{
data(geospiza)
attach(geospiza) 

# Multivariage continuous
drop.tip(geospiza.tree, "olivacea")->g.tree


s<-ic.sigma(g.tree, geospiza.data)
sims<-sim.char(g.tree, s, 100)

# Discrete character evolution
q<-list(rbind(c(-.01, .01), c(.01, -.01)))
sims<-sim.char(geospiza.tree, q, model="discrete", n=100)
  
}

\keyword{ arith }
