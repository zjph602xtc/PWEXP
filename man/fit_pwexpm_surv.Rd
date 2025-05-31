\name{pwexpm}
\alias{pwexpm}
\alias{print.pwexpm}
\alias{plot.pwexpm}
\title{Fit the Piecewise Exponential Distribution}
\description{
Fit the piecewise exponential distribution with right censoring data. User can specifity all breakpoints, some of the breakpoints or let the function estimate the breakpoints.
}
\usage{
pwexpm(Surv, data, breakpoint=NULL, nbreak=0, exclude_int=NULL, min_pt_tail=5,
        max_set=10000, seed=1818, trace=FALSE, optimizer='mle', tol=1e-4)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{Surv}{a \code{\link[survival]{Surv}} object indicating event time and status.}
  \item{data}{a data frame in which to interpret the variables named in the \code{Surv} argument.}
  \item{breakpoint}{fixed breakpoints. Pre-specifity some breakpionts. The maximum value must be earlier than the last event time.}
  \item{nbreak}{total number of breakpoints in the model. This number includes the points specified in \code{breakpoint}. If \code{nbreak=NULL}, then \code{nbreak=ceiling(8*(# unique events)^0.2)}.}
  \item{exclude_int}{an interval that excludes any estimated breakpoints (e.g., \code{exclude_int=c(10,Inf)} will exclude any estimated breakpoints after t=10). See details.}
  \item{min_pt_tail}{the minimum number of events used for estimating the tail (the hazard rate of the last piece). See details.}
  \item{max_set}{maximum estimated combination of breakpoints. }
  \item{seed}{a random seed. Do not set seed if \code{seed=NULL}.}
  \item{trace}{(internal use) logical; if TRUE, the returned data frame contains the log-likelihood of all possible breakpoints instead of the one with maximum likelihood. }
   \item{optimizer}{one of the optimizers: \code{mle}, \code{ols}, or \code{hybrid}. }
   \item{tol}{the minimum allowed gap between two breakpoints. The gap is calculated as \code{(max(time)-min(time))*tol}. Keep it as default in most cases. }
}
\details{
See webpage \url{https://zjph602xtc.github.io/PWEXP/} for a detailed description of the model and optimizers.

If user specifies \code{breakpoint}, we will check the values to make the model identifiable. Any breakpoints after the last event time will be removed;  Any breakpoints before the first event time will be removed; a mid-point will be used if there are NO events between  any two concesutive breakpoints. A warning will be given.

If user sets \code{nbreak=NULL}, then the function will automatically apply \code{nbreak=ceiling(8*(# unique events)^0.2)}. This empirical number of breakpoints is for the reference below, and it may be too large in many cases.

Argument \code{exclude_int} is a vector of two values such as \code{exclude_int=c(a, b)} (\code{b} can be \code{Inf}). It defines an interval that excludes any estimated breakpoints. It is helpful when excluding breakpoints that are too close to the tail.

In order to obtain a more robust hazard rate estimation of the tail, user can set \code{min_pt_tail} to the minimum number of events for estimating the tail (last piece of the piecewise exponential). It only works for \code{optimizer='mle'}.
}
\value{
A object of class "\code{pwexpm}" is a list containing the following components:
  \item{brk}{estimated breakpoints.}
  \item{lam}{estimated piecewise hazard rates.}
  \item{logLik}{the log-likelihood of the model.}
  \item{AIC}{the Akaike information criterion of the model.}
  \item{BIC}{the Bayesian information criterion of the model.}
  \item{para}{the parameters used to estimate the model.}
The generic accessor functions \code{AIC}, \code{BIC}, \code{logLik} can be used to extract various useful statistics from \code{pwexpm}.
The \code{plot} function can be used to make a simple plot for \code{pwexpm}.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{pwexpm_fit}}, \code{\link{boot.pwexpm}}, \code{\link{cv.pwexpm}}}

\references{
Muller, H. G., & Wang, J. L. (1994). Hazard rate estimation under random censoring with varying kernels and bandwidths. Biometrics, 61-76.
}
\examples{
event_dist <- function(n)rpwexpm(n, rate=c(0.1, 0.01, 0.2), breakpoint=c(5,14))
dat <- simdata(rand_rate = 20, total_sample = 1000, drop_rate = 0.03,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))
cut <- quantile(dat$randT, 0.8)
train <- cut_dat(var_randT = 'randT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')

fit_a0 <- pwexpm(Surv(followT, event), data=train, breakpoint = c(5,14))
fit_a1 <- pwexpm(Surv(followT, event), data=train, nbreak = 2, breakpoint = 14)
fit_b0 <- pwexpm(Surv(followT, event), data=train, nbreak = 0)
fit_b1 <- pwexpm(Surv(followT, event), data=train, nbreak = 1)
fit_b2 <- pwexpm(Surv(followT, event), data=train, nbreak = 2)
}
