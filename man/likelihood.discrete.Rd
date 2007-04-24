\name{likelihood.discrete}
\alias{likelihood.discrete}
\title{ Likelihood for discrete character evolution }
\description{
 Calculates the likelihood score for a set of discrete character data under a specified model
}
\usage{
likelihood.discrete(phy, tip.data, q, delta=1, lambda=1, endf=1, break.point=0, f=1, rtt.rescale=0, total.rescale=F)
}
\arguments{
  \item{phy}{ Phylogenetic tree }
  \item{tip.data}{ Tip data }
  \item{q}{ Instantaneous rate of change }
  \item{delta}{ Tree scaling parameter }
  \item{lambda}{ Tree scaling parameter }
  \item{endf}{ Tree scaling parameter }
  \item{break.point}{ Tree scaling parameter }
  \item{f}{ Tree scaling parameter }
  \item{rtt.rescale}{ Rescale root-to-tip distance? }
  \item{total.rescale}{ Rescale total tree length? }


}
\details{
}
\value{
 	-loglnL
}
\references{ Felsenstein pruning }
\author{ Luke J. Harmon }
\note{ 
}
\seealso{  }
\examples{
data(geospiza)
attach(geospiza)
drop.tip(geospiza.tree, "olivacea")->g.tree

gb<-as.factor(geospiza.data[,1]>4.2)
names(gb)<-rownames(geospiza.data)

likelihood.discrete(g.tree, gb, q=1)
likelihood.discrete(g.tree, gb, q=1, delta=0.5)
likelihood.discrete(g.tree, gb, q=1, lambda=0.5)
likelihood.discrete(g.tree, gb, q=1, endf=0.5)
likelihood.discrete(g.tree, gb, q=1, break.point=0.5, endf=0.2)



}

\keyword{ arith }
