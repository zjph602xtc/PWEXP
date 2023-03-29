\name{plot_survival}
\alias{plot_survival}
\alias{plot_survival.default}
\alias{plot_survival.pwexp.fit}
\title{Plot Survival Curve}
\description{
  Plot KM curve  with right censoring data or the survival curve of a fitted piecewise exponential model.
}
\usage{
\method{plot_survival}{default}(time, event, add=FALSE, conf.int=FALSE, mark.time=TRUE,
       lwd=2, xlab='Follow-up time', ylab='Survival function', ...)
\method{plot_survival}{pwexp.fit}(object, add=TRUE, show_breakpoint=TRUE,
       breakpoint_par=NULL, ...)
\method{plot_survival}{boot.pwexp.fit}(object, add=TRUE, alpah=0.05, show_breakpoint=TRUE,
       breakpoint_par=NULL, show_CI=TRUE, CI_par=NULL, ...)
}
\arguments{
 \item{time}{observed time from enrollment. For right censored data, this is the follow-up time. }
 \item{event}{the status indicator, normally 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event). }
 \item{object}{an object of class \code{pwexp.fit} or \code{boot.pwexp.fit} returned by the \code{\link{pwexp.fit}} or \code{\link{boot.pwexp.fit}} function.}
 \item{add}{logical; if TRUE add lines to current plot.}
 \item{show_breakpoint}{logical; if TRUE add vertial dashed lines to indicate breakpoints.}
 \item{breakpoint_par}{a list of parameters to control the apperance of vertical lines of breakpoionts. The values pass to \code{\link[graphics]{abline}}.}
 \item{alpha}{the significance level of the confidence interval. }
  \item{show_CI}{logical; if TRUE add confidence interval of the estimated curve. For KM esitmator, use \code{conf.int=TRUE} to show CI band.}
 \item{CI_par}{a list of parameters to control the apperance of lines of confidence intervals. The values pass to \code{\link[graphics]{lines}}.}
 \item{...}{other arguments (e.g., \code{conf.int}, \code{mark.time}, \code{lwd}, etc.) are passed over to \code{\link[survival]{plot.survfit}} (default method) or \code{\link[base]{plot}} (for class \code{pwexp.fit}).}
}
\details{
For the default method, this a  wrapper of \code{\link[survival]{plot.survfit}} function to plot right censoring data.

For class \code{pwexp.fit}, parameters in \code{...} are passed to \code{\link[base]{plot}} function to control the appearance of the survival curve; parameters in \code{breakpoint_par} are passed to \code{\link[graphics]{abline}} function to control the appearance of vertical lines of breakpoints. See examples for usage.

For class \code{boot.pwexp.fit}, parameters in \code{...} are passed to \code{\link[base]{plot}} function to control the appearance of the survival curve; parameters in \code{breakpoint_par} are passed to \code{\link[graphics]{abline}} function to control the appearance of vertical lines of breakpoints; parameters in \code{CI_par} are passed to \code{\link[graphics]{lines}} function to control the appearance of confidence intervals. See examples for usage.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}

\examples{
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(enroll_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))

# plot_survival(dat$followT, dat$event, xlim=c(0,40))

fit_res <- pwexp.fit(dat$followT, dat$event, nbreak = 2)
# plot_survival(fit_res, col='red', lwd=3, breakpoint_par = list(col='grey', lwd=2.5))
}

\seealso{\code{\link{plot_event}}}
