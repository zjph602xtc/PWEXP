\name{cv.pwexp.fit}
\alias{cv.pwexp.fit}
\alias{cv.pwexp.fit.default}
\alias{cv.pwexp.fit.pwexp.fit}
\title{Cross Validate a Piecewise Exponential Model}
\description{
Cross Validate a existing piecewise exponential model.
}
\usage{
\method{cv.pwexp.fit}{default}(time, event, nfold=5, nsim=100, breakpoint=NULL,
             nbreak=0, max_set=1000, seed=1818)
\method{cv.pwexp.fit}{pwexp.fit}(object, nfold=5, nsim=100, max_set=1000, seed=1818)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{pwexp.fit} returned by the \code{\link{pwexp.fit}} function.}
  \item{time}{observed time from enrollment. See \code{\link{pwexp.fit}}.}
  \item{event}{the status indicator.  See \code{\link{pwexp.fit}}.}
  \item{nfold}{the number of folds used in CV.}
  \item{nsim}{the number of simulations.}
  \item{breakpoint}{pre-specified breakpoints.   See \code{\link{pwexp.fit}}.}
  \item{nbreak}{total number of breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{max_set}{maximum estimated combination of breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{seed}{a random seed.}
}
\details{
Use cross validation obtain the prediction log likelihood.
}
\value{
A vector of length \code{nsim} containing the CV log likelihood in each round of simulation.
}


\examples{
\dontrun{
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(enroll_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))

cv0 <- cv.pwexp.fit(dat$followT, dat$event, nsim = 100, nbreak = 0)
cv1 <- cv.pwexp.fit(dat$followT, dat$event, nsim = 100, nbreak = 1)
cv2 <- cv.pwexp.fit(dat$followT, dat$event, nsim = 100, nbreak = 2)
sapply(list(cv0,cv1,cv2), median)
}
}

\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{pwexp.fit}}}
