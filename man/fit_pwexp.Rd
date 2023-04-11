\name{pwexp.fit}
\alias{pwexp.fit}
\title{Fit the Piecewise Exponential Distribution}
\description{
Fit the piecewise exponential distribution with right censoring data. User can specifity all breakpoints, some of the breakpoints or let the function estimate the breakpoints.
}
\usage{
pwexp.fit(time, event, breakpoint=NULL, nbreak=0, exclude_int=NULL, max_set=10000,
        seed=1818, trace=FALSE, optimizer='mle')
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{time}{observed time from randomization. For right censored data, this is the follow-up time.}
  \item{event}{the status indicator, normally 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event). }
  \item{breakpoint}{fixed breakpoints. Pre-specifity some breakpionts. The maximum value must be earlier than the last event time.}
  \item{nbreak}{total number of breakpoints in the model. This number includes the points specified in \code{breakpoint}.}
  \item{exclude_int}{an interval that excludes any estimated breakpoints (e.g., \code{exclude_int=c(10,Inf)} will exclude any estimated breakpoints after t=10) }
  \item{max_set}{maximum estimated combination of breakpoints. }
  \item{seed}{a random seed.}
  \item{trace}{(internal use) logical; if TRUE, the returned data frame contains the log-likelihood of all possible breakpoints instead of the one with maximum likelihood. }
   \item{optimizer}{one of the optimizers: \code{mle}, \code{ols}, or \code{hybrid}. }
}
\details{
See webpage \url{https://zjph602xtc.github.io/PWEXP/} for a detailed description of the model and optimizers.

If user specifies \code{breakpoint}, we will check the values to make the model identifiable. Any breakpoints after the last event time will be removed;  Any breakpoints before the first event time will be removed; a mid-point will be used if there are NO events between  any two concesutive breakpoints. A warning will be given.

Argument \code{exclude_int} is a vector of two values such as \code{exclude_int=c(a, b)} (\code{b} can be \code{Inf}). It defines an interval that excludes any estimated breakpoints. It is helpful when excluding breakpoints that are too close to the tail.
}
\value{
A data frame (\code{res}) containing these columns:
  \item{brk1, ..., brkx}{estimated breakpoints. The \code{attr(res,'brk')} can extract the vector of breakpoint from the model (\code{res} is the returned model from \code{pwexp.fit}). }
  \item{lam1, ..., lamx}{estimated piecewise hazard rates. The \code{attr(res,'lam')} can extract the vector of hazard rates from the model (\code{res} is the returned model from \code{pwexp.fit}). }
  \item{likelihood}{the log-likelihood of the model.}
  \item{AIC}{the Akaike information criterion of the model.}
  \item{BIC}{the Bayesian information criterion of the model.}
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{boot.pwexp.fit}}, \code{\link{cv.pwexp.fit}}}

\examples{
event_dist <- function(n)rpwexp(n, rate=c(0.1, 0.01, 0.2), breakpoint=c(5,14))
dat <- simdata(rand_rate = 20, total_sample = 1000, drop_rate = 0.03,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))
cut <- quantile(dat$randT, 0.8)
train <- cut_dat(var_randT = 'randT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')

fit_a0 <- pwexp.fit(train$followT, train$event, breakpoint = c(5,14))
fit_a1 <- pwexp.fit(train$followT, train$event, nbreak = 2, breakpoint = c(14))
fit_b0 <- pwexp.fit(train$followT, train$event, nbreak = 0)
fit_b1 <- pwexp.fit(train$followT, train$event, nbreak = 1)
fit_b2 <- pwexp.fit(train$followT, train$event, nbreak = 2)
}
