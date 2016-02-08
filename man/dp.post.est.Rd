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
\details{Performs probability mass function estimation under th following model
\deqn{y_i \mid P \sim  P, i=1, \dots, n}{y_i | P ~ P, i=1, ..., n}
\deqn{P \sim DP(\alpha, P_0),}{P ~ Pi(\alpha, P_0),}
where \eqn{P_0} is Poisson with mean \code{lambda}.}
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