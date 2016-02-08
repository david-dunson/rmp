\name{npmp}
\alias{dpmpoiss}
\alias{npmp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DP mixtures of Poissons
}
\description{
Performs probability mass function estimation under nonparametric mixture of Poisson kernels. }
\usage{
npmp(y, k, nrep, nb, alpha=1, theta=alpha, sigma=0, 
mixing_hyperprior= FALSE, basemeasure_hyperprior = FALSE, mixing_type="DP", 
algo="slice",  prior="gamma", a, b, a_a=1, b_a=1, lb=NULL, ub=NULL, 
print = 1, ndisplay = nrep/4, plot.it = FALSE, pdfwrite = FALSE, ... )

dpmpoiss(y, k, nrep, nb, alpha = 1, a, b, lb = NULL, ub = NULL, print = 1,
ndisplay = nrep/4, plot.it = FALSE, pdfwrite = FALSE, ...)
}
\arguments{
  \item{y}{
Vector of count data
}
  \item{k}{
Truncation level for the number of cluster in the mixture. Default is \code{length(ydis)}.
}
  \item{nrep}{
Number of MCMC iterations
}
  \item{nb}{
Number of burn-in iteration in the MCMC to discard
}
  \item{alpha}{
Value of the precision parameter of the Dirichlet process prior
}
  \item{theta}{
Value of the strength parameter of the Two-parameters-Poisson-Dirichlet process prior
}
  \item{sigma}{
Value of the discount parameter of the Two-parameters-Poisson-Dirichlet process prior
}
  \item{mixing_hyperprior}{
Logical. If TRUE \code{alpha} is random with gamma hyperprior
}
  \item{basemeasure_hyperprior}{
Logical. If TRUE also the parameters of the base measure are random, see details below.
}
  \item{mixing_type}{
Type of mixing distribution. Default is "DP" for Dirichlet process but also "2PD" for Two-parameters-Poisson-Dirichlet process is allowed.
}
 \item{algo}{
Type of algorithm. Current choices are: slice sampler (\code{algo="slice"}) or polya-urn-type sampler (\code{algo="polya-urn"}).
}
  \item{prior}{
String for the base measure prior. Default is "gamma" for lambda ~ Gamma(a,b). The other choice is "normal" for exp(lambda) ~ N(a,b)}
  \item{a}{
Shape (mean) hyperparameter for the gamma (normal) base measure 
}
  \item{b}{
Scale (sd) hyperparameter for the gamma (normal) prior
}
  \item{a_a}{
Shape hyperparameter for \code{alpha}
}
  \item{b_a}{
Scale  hyperparameter for \code{alpha}
}
  \item{lb}{
Scalar integer. Lower bound for the argument of the pmf. Default is \code{max(0,min(ydis)-10)}.
}
  \item{ub}{
Scalar integer. Upper bound for the argument of the pmf. Default is \code{max(ydis)+10}.
}
  \item{print}{
Vector of integers (from 1 to 5) indicating whether to print each step of the Gibbs sampler.  Specifically, 1 for current iteration, 2 for the DP cluster allocation, 3 for the posterior parameters of the mixture components, 4 for the precision of the DP, 5 for the posterior pmf.
}
  \item{ndisplay}{
Scalar integer. It gives the number of iterations to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out)
}
  \item{plot.it}{
Logical, default FALSE. If TRUE a plot with empirical and estimated posterior probability mass functions is plotted.
}
  \item{pdfwrite}{
Logical, default FALSE. If TRUE a pdf file is written in the current working directory. Traceplots and other posterior quantities are drown.
}
  \item{...}{
Additional arguments (for future implemetantions).
}
}
\details{The function \code{npmp} performs probability mass function estimation under nonparametric mixture of Poisson kernels, i.e.
\deqn{y_i \mid \lambda_i \sim  \mbox{Poi}(\lambda_i), i=1, \dots, n}{y_i | \lambda_i ~ Poi(\lambda_i), i=1, ..., n}
\deqn{\lambda_i \mid G \sim G}{\lambda_i | G ~ G}
\deqn{G \sim \Pi(P_0),}{G ~ Pi(P_0),}
where \eqn{\Pi} is a nonparametric prior (Dirichlet process or Two-parameters-Poisson-Dirichlet process) with base measure \eqn{P_0}. The function \code{dppoiss} is a wrapper to \code{npmp} with \code{mixing_type="DP"} for back portability with version 1.0 of the package. The main part of the code is written in C language to gain computational speed. Plots and posterior summaries are in plain R code. From version 2.0 on, the blocked Gibbs sampler has been removed in place of slice samper (Kalli et al., 2011) and polya-urn sampler. Two different base measures \eqn{P_0} are implemented: \code{prior="gamma"} and \code{prior="normal"} for 
	\deqn{
\lambda_h \sim \mbox{Gamma}(a,b), \log(\lambda_h) \sim N(a,b), h=1,\dots
}{%
\lambda_h ~ Gamma(a,b), log(\lambda_h) ~ N(a,b), h=1,...
}
 respectively.}
\value{
\item{name}{Name of the model}
\item{mixing_type}{Name of the mixing prior}
\item{mcmc}{Quantities about MCMC sampling}
\item{mcmc.chains}{MCMC chains of the parameters}
\item{pmf}{A list containing several quantities related to the probability mass function (emprical pmf, posterior mean pmf and pointwise 95\% credible intervals) computed for the values from \code{lb} to \code{ub}}
\item{parameters}{A list containing the posterior mean of the cluster specific parameters (be careful of label-switching problems)}
\item{clustering}{A list containing posterior quantities related to the clustering structure of the data} 
}
\references{
Canale, A. and Dunson, D. B. (2011), "Bayesian Kernel Mixtures for Counts", \emph{Journal of American Statistical Association}, \bold{106}, 1528-1539.

Kalli, M., Griffin, J., and Walker, S. (2011), "Slice sampling mixture models," Statistics and Computing, \bold{21}, 93-105.
}
\author{
R code and porting by A. Canale, C code by A. Canale with minor contributions by N. Lunardon. 
}
\seealso{
rmg
}
\examples{
\dontrun{
data(ethylene)
y <- tapply(ethylene$impl,FUN=mean,INDEX=ethylene$id)
z <- tapply(ethylene$dose,FUN=mean,INDEX=ethylene$id)

# Estimate the pmf of the number of implants in the control group
y0  <- y[z==0]
pmf.control = dpmpoiss(y0, k=20, nrep=11000, nb=1000, alpha=1, a=1, b=1,
lb=5, ub=24, plot.it=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
