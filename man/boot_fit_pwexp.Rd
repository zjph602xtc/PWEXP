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
               exclude_int=NULL, max_set=1000, seed=1818, optimizer='mle')
\method{boot.pwexp.fit}{pwexp.fit}(object, nsim=100, max_set=1000, seed=1818, optimizer='mle')
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{pwexp.fit} returned by the \code{\link{pwexp.fit}} function.}
  \item{time}{observed time from randomization. See \code{\link{pwexp.fit}}.}
  \item{event}{the status indicator.  See \code{\link{pwexp.fit}}.}
  \item{nsim}{the number of repeated bootstraping.}
  \item{breakpoint}{pre-specified breakpoints.   See \code{\link{pwexp.fit}}.}
  \item{nbreak}{total number of breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{exclude_int}{an interval that excludes any estimated breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{max_set}{maximum estimated combination of breakpoints.  See \code{\link{pwexp.fit}}.}
  \item{seed}{a random seed.}
  \item{optimizer}{one of the optimizers: \code{mle}, \code{ols}, or \code{hybrid}. See \code{\link{pwexp.fit}}.}
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
\dontrun{
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(rand_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))
cut <- quantile(dat$randT, 0.8)
train <- cut_dat(var_randT = 'randT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')

fit_res3 <- pwexp.fit(train$followT, train$event, nbreak = 2)
fit_res_boot <- boot.pwexp.fit(fit_res3, nsim = 100)
plot_survival(train$followT, train$event, xlim=c(0,40))
plot_survival(fit_res_boot, col='red', CI_par = list(col='red'))
brk_ci <- apply(attr(fit_res_boot, 'brk'), 2, function(x)quantile(x,c(0.025,0.975)))
abline(v=brk_ci, col='grey', lwd=2)
}
}
