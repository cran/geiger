\name{fit.discrete}
\alias{fit.discrete}
\title{ Model fitting for discrete comparative data }
\description{
	Fits macroevolutionary models to phylogenetic trees
}
\usage{
fit.discrete(phy, data, data.names=NULL, lambda=FALSE, delta=FALSE, kappa=FALSE, linearchange=FALSE, exponentialchange=FALSE, tworate=FALSE, start.rate=1.0)
}
\arguments{
  \item{phy}{ object of type phylo }
  \item{data}{ Data vector (one trait) or matrix (multiple traits) }
  \item{data.names}{ Tip names for data vector that match tree species; ignored if data includes names}
  \item{lambda}{ Estimate lambda? }
  \item{delta}{ Estimate delta? }
  \item{kappa}{ Estimate kappa? }
  \item{linearchange}{ Fit model where rate changes linearly through time? }
  \item{exponentialchange}{ Fit model where rate changes exponential through time? }
  \item{tworate}{ Fit model where rate changes at some time point in the tree? }
  \item{start.rate}{ Starting guess for parameters }
}
\details{
  This function fits various likelihood models for discrete character evolution.   Likelihood is maximized 
  using the r function nlm. All of the models are 
  continuous-time Markov models of trait evolution (see Yang 2006 for a good general discussion of this 
  type of model).  The function can handle traits with any number of character states, although it has a 
  major limitation: it can only implement a model where all transitions are equally likely.  That is, the 
  approach assumes that traits are evolving in a Jukes-Cantor like model, where changes from any character 
  state to any other occur at the same rate, q.  There are two other parameters that represent models where 
  transition rates vary across the tree: 
 \item{lambda}{ Pagel's lambda; multiplies all internal branches of the tree by lambda, leaving tip branches as their original length.}
 \item{kappa}{ Pagel's kappa; raises all branch lengths to the power kappa.  As kappa approaches zero, the model becomes speciational.}
 \item{delta}{ Pagel's delta; raises all node depths to the power delta. If delta is less than one, evolution in concentrated early in the tree; delta > 1 concentrates evolution towards the tips.}
 \item{linearchange}{ A model where the rate of evolution changes linearly through time.
						The rate at the present day is supplied as endRate, and only needs to be specified as a ratio relative to the initial rate
						at the root of the tree.  Rates change linearly: r(t) = ro + a * t, where ro is the inital rate and a is the slope determined using endRate. 
						If endRate is greater than 1, evolution gradually speeds up; if endrate is less than 1, evolution gradually slows down. If endRate = 1, the model is a constant-
						rate Brownian motion model (so the tree is returned unchanged).}
 \item{exponentialchange}{ A model where the rate of evolution changes exponentially through time.
						The rate at the present day is supplied as endRate, in exactly the same way as in linearchange.tree.
						The difference for this function is that rates grow or decay exponentially rather than linearly.
						The equation is r(t) = ro * exp(a * t), where ro is the inital rate and a is determined using endRate.}
 \item{tworate}{ A model that effectively changes the rate of evolution at some point in time, specified by breakPoint.
						The rate after the breakpoint is supplied as endRate, and only needs to be specified as a ratio relative to the initial rate.  
						Thus, if endRate is greater than 1, evolution speeds up, all branches towards the tips are made longer,  
						if endrate is less than 1, branches towards the tips are made shorter. If endRate = 1, the model is a constant-
						rate model.}
 \item{start.rate}{ Starting values for the parameters you're estimating.  For a basic model (all parameters FALSE), p should contain one
  						number; for other models, more starting parameters are needed.}
}
\value{
 	Returns maximum likelihood value for q and selected parameters, along with the likelihood score.
	Sometimes, parameters are confounded!
}
\references{ 
Yang, Z. 2006. Computational Molecular Evolution. Oxford University Press: Oxford.
}
\author{ Luke J. Harmon and R. E. Glor }
\note{ 
}
\seealso{ }
\examples{
data(geospiza)
attach(geospiza)

gb<-as.factor(geospiza.data[,1]>4.2)
names(gb)<-rownames(geospiza.data)

fit.discrete(geospiza.tree, gb)
fit.discrete(geospiza.tree, gb, lambda=TRUE)
fit.discrete(geospiza.tree, gb, delta=TRUE)
fit.discrete(geospiza.tree, gb, kappa=TRUE)

fit.discrete(geospiza.tree, gb, linearchange=TRUE)
fit.discrete(geospiza.tree, gb, exponentialchange=TRUE)
fit.discrete(geospiza.tree, gb, tworate=TRUE)




}
\keyword{ arith }
