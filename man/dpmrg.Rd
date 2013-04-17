\name{dpmrg}
\alias{dpmrg}
\title{
DP mixture of rounded Gaussians
}
\description{
Performs Bayesian probability mass function estimation under rounded mixture of Gaussian prior.
}
\usage{
dpmrg(ydis, k, nrep, nb, alpha, alpha_r = FALSE, mu0 = mean(ydis), 
kap = var(ydis), atau, btau, a_a = 1, b_a = 1,  lb = NULL, ub = NULL, print = 1, 
ndisplay = nrep/4, plot = FALSE, pdfwrite = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ydis}{
Vector of count data
}
  \item{k}{
Truncation level for the number of cluster in the DP mixture
}
  \item{nrep}{
Number of MCMC iterations
}
  \item{nb}{
Number of burn-in iteration in the MCMC to discard
}
  \item{alpha}{
Value of the precision parameter of the DP
}
  \item{alpha_r}{
Logical. If TRUE \code{alpha} is also random and assigned a hyoerprior
}
  \item{mu0}{
Location hyperparameter for the latent rounded Gaussian base measure
}
  \item{kap}{
Precision hyperparameter for the latent rounded Gaussian base measure
}
  \item{atau}{
Shape hyperparameter for the Gamma distribution
}
  \item{btau}{
Scale hyperparameter for the Gamma distribution 
}
  \item{a_a}{
Shape hyperparameter for the Gamma distribution for \code{alpha} 
}
  \item{b_a}{
Scale hyperparameter for the Gamma distribution for \code{alpha} 
}
  \item{lb}{
Scalar integer. Lower bound for the argument of the pmf. Default is \code{max(0,min(ydis)-10)}.
}
  \item{ub}{
Scalar integer. Upper bound for the argument of the pmf. Default is \code{max(ydis)+10}.
}
  \item{print}{
Vector of integers (from 1 to 6) indicating whether to print each step of the Gibbs sampler.  Specifically, 1 for current iteration, 2 for the data augmentation step simulating the latent continuous variables, 3 for the DP cluster allocation, 4 for the posterior parameters of the mixture components, 5 for the precision of the DP, 6 for the posterior pmf.
}
  \item{ndisplay}{
Scalar integer. It gives the number of iterations to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out)
}
  \item{plot}{
Logical, default FALSE. If TRUE a plot with empirical and estimated posterior probability mass functions is plotted.
}
  \item{pdfwrite}{
Logical, default FALSE. If TRUE a pdf file is written in the current working directory. Traceplots and other posterior quantities are drown.
}
  \item{...}{
Additional arguments (for future implemetantions).
}
}
\details{
The function performs Bayesian probability mass function estimation under the mixture model of Canale and Dunson (2011).
The main part of the code is written in C language to gain computational speed. Plots and posterior summaries are in plain R code.
}
\value{
\item{name}{The name of the function}
\item{mcmc}{Quantities about MCMC sampling}
\item{pmf}{A list containing several quantities related to the probability mass function (emprical pmf, posterior mean pmf and pointwise 95\% credible intervals) computed for the values from \code{lb} to \code{ub}}
\item{parameters}{A list containing the posterior mean of the parameters}
\item{clustering}{A list containing posterior quantities related to the clustering structure of the data} 
}
\references{
Canale, A. and Dunson, D. B. (2011), "Bayesian Kernel Mixtures for Counts", \emph{Journal of American Statistical Association}, \bold{106}, 1528-1539.
}
\author{
C code by A. Canale and N. Lunardon, R code and porting by A. Canale
}
\seealso{
dpmpoiss
}
\examples{
\dontrun{
data(ethylene)
y <- tapply(ethylene$impl,FUN=mean,INDEX=ethylene$id)
z <- tapply(ethylene$dose,FUN=mean,INDEX=ethylene$id)

# Estimate the pmf of the number of implants in the control group
y0  <- y[z==0]
pmf.control = dpmrg(y0, k=20, nrep=11000, nb=1000, alpha=1, atau=1, btau=1,
lb=5, ub=21, plot= TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
