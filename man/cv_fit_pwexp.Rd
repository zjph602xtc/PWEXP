\name{cv.pwexp.fit}
\alias{cv.pwexp.fit.default}
\alias{cv.pwexp.fit.pwexp.fit}
\alias{}
\title{Cross Validate a Piecewise Exponential Model}
\description{
BCross Validate a existing piecewise exponential model.
}
\usage{
\method{cv.pwexp.fit}{default}(time, event, nsim=100, breakpoint=NULL, nbreak=0,
               max_set=1000, seed=1818)
\method{cv.pwexp.fit}{pwexp.fit}(object, nsim=100, max_set=1000, seed=1818)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{pwexp.fit} returned by the \code{\link{pwexp.fit}} function.}
  \item{time}{observed time from enrollment. See \code{\link{pwexp.fit}}.}
  \item{event}{the status indicator.  See \code{\link{pwexp.fit}}.}
  \item{breakpoint}{pre-specified breakpoints.   See \code{\link{pwexp.fit}}.}
  \item{nbreak}{total number of breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{max_set}{maximum estimated combination of breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{seed}{a random seed.}
}
\details{
Use bootstrap to repeatdly call \code{\link{pwexp.fit}} to estimate the uncertainty of parameters.
}
\value{
A data frame containing these columns:
  \item{brk1, ..., brkx}{estimated breakpoints. The \code{attr(res,'brk')} can extract the result of breakpoint from each bootstrap sample (\code{res} is the returned model from \code{boot.pwexp.fit}). }
  \item{lam1, ..., lamx}{estimated piecewise hazard rates. The \code{attr(res,'lam')} can extract the result of hazard rates from bootstrap sample (\code{res} is the returned model from \code{boot.pwexp.fit}). }
  \item{likelihood}{the log-likelihood of the model.}
  \item{AIC}{the Akaike information criterion of the model.}
  \item{BIC}{the Bayesian information criterion of the model.}
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{boot.pwexp.fit}}, \code{\link{cv.pwexp.fit}}}

\examples{

}
