\name{dp.post.est}
\alias{dp.post.est}
\title{
Posterior probability mass function estimation with DP prior
}
\description{
Performs Bayesian probability mass function estimation under DP prior with Poisson base measure.
}
\usage{
dp.post.est(x, y, alpha, lambda)
}
\arguments{
  \item{x}{
Values on which to compute the pmf.
}
  \item{y}{
Vector of observed data.
}
  \item{alpha}{
DP precision parameter
}
  \item{lambda}{
Mean parameter for the Poisson base measure.
}
}

\value{
A vector of size \code{length(x)} containing the probability masses 
}
\references{
Carota, C., and Parmigiani, G. (2002), ``Semiparametric Regression for Count Data,'' \emph{Biometrika}, \bold{89}, 265--281.
}
\author{
Antonio Canale
}

\examples{
data(ethylene)
y <- tapply(ethylene$impl,FUN=mean,INDEX=ethylene$id)
z <- tapply(ethylene$dose,FUN=mean,INDEX=ethylene$id)

# Estimate the pmf of the number of implants in the control group
y0  <- y[z==0]
pmf.control = dp.post.est(0:30, y0, alpha = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
