\name{phy.anova}
\alias{phy.anova}
\alias{phy.manova}
\title{Phylogenetic ANOVA and MANOVA}
\description{
	Calculates ANOVA or MANOVA, and returns p-value based on Brownian motion simulations
}
\usage{
phy.anova(y, factor, phy, nsim=1000)
phy.manova(y, factor, phy, nsim=1000, test="Wilks")

}
\arguments{
  \item{y}{ Independent variable(s) - continuous }
  \item{factor}{ Dependent variable - factor }
  \item{phy}{ Phylogenetic tree }
  \item{nsim}{ Number of simulations for calculating p-value }
  \item{test}{ The name of the test statistic to be used. Partial matching is used so the name can be abbreviated. }
}
\details{
Simulations are run under a Brownian motion model.  For MANOVA the estimates vcv matrix from 
ic.sigma is used as a model.

You can specify the test statistic for the summary table. Wilks' statistic is most popular in the literature.
}
\value{
	Standard ANOVA or MANOVA table and p-value based on simulations
	}
\references{ Garland }
\author{ Luke J. Harmon }
\note{ 

}
\seealso{  }
\examples{

data(geospiza)
attach(geospiza)

drop.tip(geospiza.tree, "olivacea")->g.tree

f<-as.factor(c(rep(0, 7), rep(1, 6)))


phy.manova(geospiza.data, f, g.tree)

x<-geospiza.data[,1]
names(x)<-rownames(geospiza.data)


phy.anova(x, f, g.tree) 

}

\keyword{ arith }
