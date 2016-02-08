\name{ethylene}
\alias{ethylene}
\docType{data}
\title{
  Developmental toxicity study of ethylene glycol in mice
}
\description{
Data from the developmental toxicity study of ethylene glycol in mice conducted by the National Toxicology Program. Pregnant mice were assigned to dose groups of 0, 750, 1,500 or 3,000 mg/kg per day, with the number of implants measured for each mouse at the end of the experiment. Group sizes are 25, 24, 23 and 23, respectively. 
}
\usage{data(ethylene)}
\format{
  A data frame with 1192 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{identifyer for pregnant mouse}
    \item{\code{dose}}{dose of ethylene glycol}
    \item{\code{weight}}{weight of the fetus}
    \item{\code{sex}}{sex of the fetus}
    \item{\code{impl}}{number of implants in the pregnant mouce}
    \item{\code{litsz}}{size of the relative litter}
    \item{\code{malf}}{presence of malformation in the fetus}
  }
}
\references{
  Price, C. J., Kimmel, C. A., Tyl, R. W., and Marr, M. C.  (1985) "The developemental toxicity of ethylene glycol in rats and mice" \emph{Toxicological and Applied Pharmacology} \bold{81}, 113-127.
}
\examples{
data(ethylene)
implants <- tapply(ethylene$impl, FUN=mean, INDEX=ethylene$id)

summary(implants)
m <- mean(implants)
v <- var(implants)
hist(implants, main=paste("Histogram of the number of implants (Mean = ",
round(m,2), ", Var = ", round(v,2),")"))
}
\keyword{datasets}
