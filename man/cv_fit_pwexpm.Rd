\name{cv.pwexpm_fit}
\alias{cv.pwexpm_fit}
\title{Cross Validate a Piecewise Exponential Model}
\description{
Build and cross validate a piecewise exponential model.
}
\usage{
cv.pwexpm_fit(time, event, nfold=5, nsim=100, breakpoint=NULL,
             nbreak=0, exclude_int=NULL, min_pt_tail=5, max_set=1000, seed=1818,
             optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{time}{observed time from randomization.}
  \item{event}{the status indicator.  See \code{\link{pwexpm_fit}}.}
  \item{nfold}{the number of folds used in CV.}
  \item{nsim}{the number of simulations.}
  \item{breakpoint}{pre-specified breakpoints. See \code{\link{pwexpm_fit}}.}
  \item{nbreak}{total number of breakpoints.  See \code{\link{pwexpm_fit}}.}
  \item{exclude_int}{an interval that excludes any estimated breakpoints.  See \code{\link{pwexpm_fit}}.}
  \item{min_pt_tail}{the minimum number of events used for estimating the tail (the hazard rate of the last piece).  See \code{\link{pwexpm_fit}}.}
  \item{max_set}{maximum estimated combination of breakpoints.  See \code{\link{pwexpm_fit}}.}
  \item{seed}{a random seed. Do not set seed if \code{seed=NULL}.}
  \item{optimizer}{one of the optimizers: \code{mle}, \code{ols}, or \code{hybrid}. See \code{\link{pwexpm_fit}}.}
  \item{tol}{the minimum allowed gap between two breakpoints. The gap is calculated as \code{(max(time)-min(time))*tol}. Keep it as default in most cases. }
  \item{parallel}{logical. If \code{TRUE}, use \pkg{doSNOW} package to run in parallel.}
  \item{mc.core}{number of processes allowed to be run in parallel.}
  \item{...}{internal function reserved. }
}
\details{
Use cross validation obtain the prediction log likelihood.
}
\value{
A object of class "\code{cv.pwexpm}" is a numeric vector containing the CV log likelihood in each round of simulation. The \code{plot} function can be used to make a boxplot of the CV log likelihoods from \code{pwexpm}.
}


\examples{
event_dist <- function(n)rpwexpm(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(rand_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))

# here nsim=10 is for demo purpose, pls increase it in practice!!
\donttest{
cv0 <- cv.pwexpm_fit(dat$followT, dat$event, nsim = 10, nbreak = 0)
cv1 <- cv.pwexpm_fit(dat$followT, dat$event, nsim = 10, nbreak = 1)
cv2 <- cv.pwexpm_fit(dat$followT, dat$event, nsim = 10, nbreak = 2)
sapply(list(cv0,cv1,cv2), median)
}
}

\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{cv.pwexpm}}}
