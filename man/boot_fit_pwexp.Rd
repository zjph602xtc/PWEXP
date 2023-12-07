\name{boot.pwexp.fit}
\alias{boot.pwexp.fit}
\alias{boot.pwexp.fit.default}
\alias{boot.pwexp.fit.pwexp.fit}
\title{Bootstrap a Piecewise Exponential Model}
\description{
Bootstrap a existing piecewise exponential model or build a piecewise exponential model with bootstrapping.
}
\usage{
\method{boot.pwexp.fit}{default}(time, event, nsim=100, breakpoint=NULL, nbreak=0,
               exclude_int=NULL, min_pt_tail=5, max_set=1000, seed=1818,
               optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...)
\method{boot.pwexp.fit}{pwexp.fit}(time, nsim=100, max_set=1000, seed=1818,
               optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{time}{observed time from randomization or a \code{\link{pwexp.fit}} object.}
  \item{event}{the status indicator. See \code{\link{pwexp.fit}}.}
  \item{nsim}{the number of repeated bootstraping.}
  \item{breakpoint}{pre-specified breakpoints. See \code{\link{pwexp.fit}}.}
  \item{nbreak}{total number of breakpoints. See \code{\link{pwexp.fit}}.}
  \item{exclude_int}{an interval that excludes any estimated breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{min_pt_tail}{the minimum number of events used for estimating the tail (the hazard rate of the last piece). See \code{\link{pwexp.fit}}.}
  \item{max_set}{maximum estimated combination of breakpoints. See \code{\link{pwexp.fit}}.}
  \item{seed}{a random seed.}
  \item{optimizer}{one of the optimizers: \code{mle}, \code{ols}, or \code{hybrid}. See \code{\link{pwexp.fit}}.}
  \item{tol}{the minimum allowed gap between two breakpoints. The gap is calculated as \code{(max(time)-min(time))*tol}. Keep it as default in most cases. }
  \item{parallel}{logical. If \code{TRUE}, use \pkg{doSNOW} package to run in parallel.}
  \item{mc.core}{number of processes allowed to be run in parallel.}
  \item{...}{internal function reserved. }
}
\details{
Use bootstrap to repeatdly call \code{\link{pwexp.fit}} to estimate the uncertainty of parameters.
}
\value{
A data frame (\code{res}) containing these columns:
  \item{brk1, ..., brkx}{estimated breakpoints. The \code{attr(res,'brk')} can extract the result of breakpoint from each bootstrap sample (\code{res} is the returned model from \code{boot.pwexp.fit}). }
  \item{lam1, ..., lamx}{estimated piecewise hazard rates. The \code{attr(res,'lam')} can extract the result of hazard rates from bootstrap sample (\code{res} is the returned model from \code{boot.pwexp.fit}). }
  \item{likelihood}{the log-likelihood of the model.}
  \item{AIC}{the Akaike information criterion of the model.}
  \item{BIC}{the Bayesian information criterion of the model.}
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{pwexp.fit}}}

\examples{
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(rand_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))

fit_res3 <- pwexp.fit(dat$followT, dat$event, nbreak = 2)
fit_res_boot <- boot.pwexp.fit(fit_res3, nsim = 10) # here nsim=10 is for demo purpose,
                                                    # pls increase it in practice
plot_survival(dat$followT, dat$event, xlim=c(0,40))
plot_survival(fit_res_boot, col='red', CI_par = list(col='red'))
brk_ci <- apply(attr(fit_res_boot, 'brk'), 2, function(x)quantile(x,c(0.025,0.975)))
abline(v=brk_ci, col='grey', lwd=2)
}
