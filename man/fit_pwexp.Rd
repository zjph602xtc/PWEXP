\name{pwexp.fit}
\title{Fit the Piecewise Exponential Distribution}
\description{
Fit the piecewise exponential distribution with right censoring data. User can specifity all breakpoints, some of the breakpoints and let the function estimate the breakpoints.
}
\usage{
pwexp.fit(time, event, breakpoint=NULL, nbreak=0, max_set=10000,
        seed=1818, trace=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{time}{observed time from enrollment. For right censored data, this is the follow-up time.}
  \item{event}{the status indicator, normally 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event). }
  \item{breakpoint}{Fixed breakpoints. Pre-specifity some breakpionts. The maximum value must be earlier than the last event time.}
  \item{nbreak}{total number of breakpoints in the model. This number includes the points specified in \code{breakpoint}.}
  \item{max_set}{maximum estimated combination of breakpoints. }
  \item{seed}{a random seed.}
  \item{trace}{logical; if TRUE, the returned data frame contains the log-likelihood of all possible breakpoints instead of the one with maximum likelihood. }
}
\details{
If user specifies \code{breakpoint}, we will check the values to make the model identifiable. Any breakpoints after the last event time will be removed; a mid-point will be used if there are NO events between  any two concesutive breakpoints. A warning will be given.
}
\value{
A data frame containing these columns:
  \item{brk1, ..., brkx}{estimated breakpoints. The \code{attr(res,'brk')} can extract the vector of breakpoint from the model (\code{res} is the returned model from \code{pwexp.fit}). }
  \item{lam1, ..., lamx}{estimated piecewise hazard rates. The \code{attr(res,'lam')} can extract the vector of hazard rates from the model (\code{res} is the returned model from \code{pwexp.fit}). }
  \item{likelihood} {the log-likelihood of the model.}
  \item{AIC}{the Akaike information criterion of the model.}
  \item{BIC}{the Bayesian information criterion of the model.}
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{XXX}}}

\examples{

}
