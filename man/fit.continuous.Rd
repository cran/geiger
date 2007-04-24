\name{fit.continuous}
\alias{fit.continuous}
\title{ Model fitting for continuous data }
\description{
	Fits macroevolutionary models to phylogenetic trees
}
\usage{
fit.continuous(tips, phy, lambda = FALSE, kappa = FALSE, delta = FALSE, alpha=FALSE, r = FALSE, bounds=NULL, print=TRUE, meserr=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tips}{ Tip data }
  \item{phy}{ Phylogenetic tree }
  \item{lambda}{ Estimate lambda? }
  \item{kappa}{ Estimate kappa? }
  \item{delta}{ Estimate delta? }
  \item{alpha}{ Estimate alpha? }
  \item{r}{ Estimate r? }
  \item{bounds}{ Range to constrain estimates? }
  \item{print}{ Verbose mode? }
  \item{meserr}{ Measurement error for each tip species }
}
\details{
}
\value{
	Prints out results - I think?
}
\references{ PAGEL }
\author{ Wendell Challenger and Luke J. Harmon}
\note{ 

}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)


V1 			<- log(geospiza.data[,1])
names(V1) 	<- rownames(geospiza.data)

drop.tip(geospiza.tree, "olivacea")->g.tree


#---- PRINT RESULTS
fit.continuous(V1, g.tree, print=TRUE)
 
#---- STORE RESULTS 
brown.fit <-  fit.continuous(V1, g.tree, print=FALSE)
#----------------------------------------------------
#   PHYLOGENETIC SIGNAL: FIT LAMBDA ON OFF-DAIGS
#----------------------------------------------------
fit.continuous(V1, g.tree, lambda=TRUE, print=TRUE)
#----------------------------------------------------
#    FIT TIME PROPORTIONALITY: DELTA 
#---------------------------------------------------
fit.continuous(V1, g.tree,  delta=TRUE, print=TRUE, bounds=list(delta=c(0.001, 5)))
#----------------------------------------------------
#    FIT TIME PROPORTIONALITY: KAPPA 
#---------------------------------------------------
fit.continuous(V1, g.tree,  kappa=TRUE, print=TRUE)

}

\keyword{ arith }
