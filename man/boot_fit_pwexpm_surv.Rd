\name{boot.pwexpm}
\alias{boot.pwexpm}
\alias{boot.pwexpm.default}
\alias{boot.pwexpm.pwexpm}
\alias{summary.boot.pwexpm}
\title{Bootstrap a PWE Model with Formula or Existing Model}
\description{
Bootstrap an existing piecewise exponential model or build a piecewise exponential model with bootstrapping.
}
\usage{
\method{boot.pwexpm}{default}(Surv, data, nsim=100, breakpoint=NULL, nbreak=0,
               exclude_int=NULL, min_pt_tail=5, max_set=1000, seed=1818,
               optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...)
\method{boot.pwexpm}{pwexpm}(Surv, nsim=100, max_set=1000, seed=1818,
               optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{Surv}{a \code{\link[survival]{Surv}} object indicating event time and status or a \code{\link{pwexpm}} object.}
  \item{data}{a data frame in which to interpret the variables named in the \code{Surv} argument.}
  \item{nsim}{the number of repeated bootstrapping.}
  \item{breakpoint}{pre-specified breakpoints. See \code{\link{pwexpm}}.}
  \item{nbreak}{total number of breakpoints. See \code{\link{pwexpm}}.}
  \item{exclude_int}{an interval that excludes any estimated breakpoints.  See \code{\link{pwexpm}}.}
  \item{min_pt_tail}{the minimum number of events used for estimating the tail (the hazard rate of the last piece). See \code{\link{pwexpm}}.}
  \item{max_set}{maximum estimated combination of breakpoints. See \code{\link{pwexpm}}.}
  \item{seed}{a random seed. Do not set seed if \code{seed=NULL}.}
  \item{optimizer}{one of the optimizers: \code{mle}, \code{ols}, or \code{hybrid}. See \code{\link{pwexpm}}.}
  \item{tol}{the minimum allowed gap between two breakpoints. The gap is calculated as \code{(max(time)-min(time))*tol}. Keep it as default in most cases. }
  \item{parallel}{logical. If \code{TRUE}, use \pkg{doSNOW} package to run in parallel.}
  \item{mc.core}{number of processes allowed to be run in parallel.}
  \item{...}{internal function reserved. }
}
\details{
Use bootstrap to repeatdly call \code{\link{pwexpm}} to estimate the uncertainty of parameters.
}
\value{
A object of class "\code{boot.pwexpm}" is a list containing the following components:
  \item{brk}{estimated breakpoints in each row.}
  \item{lam}{estimated piecewise hazard rates in each row.}
  \item{logLik}{the log-likelihood of the original model.}
  \item{AIC}{the Akaike information criterion of the original model.}
  \item{BIC}{the Bayesian information criterion of the original model.}
  \item{para}{the parameters used to estimate the model.}
The \code{plot} function can be used to make a simple plot for \code{boot.pwexpm}.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{boot.pwexpm_fit}}}

\examples{
event_dist <- function(n)rpwexpm(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(rand_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))

fit_res3 <- pwexpm(Surv(followT, event), data = dat, nbreak = 2)
fit_res_boot <- boot.pwexpm(fit_res3, nsim = 10) # here nsim=10 is for demo purpose,
                                                  # pls increase it in practice
plot_survival(dat$followT, dat$event, xlim=c(0,40))
plot_survival(fit_res_boot, col='red', CI_par = list(col='red'))
brk_ci <- confint(fit_res_boot, level = 0.95, parm = 'brk')
abline(v=brk_ci, col='grey', lwd=2)
}
