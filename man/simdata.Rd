\name{simdata}
\alias{simdata}
\title{Simulate Survival Data}
\description{
\code{simdata} is used to simulate a clinical trial data with time-to-event endpoints.
}
\usage{
simdata(group="Group 1", strata="Strata 1", allocation=1,
    event_lambda=NA, drop_rate=NA, death_lambda=NA, n_rand=NULL,
    rand_rate=NULL, total_sample=NULL, add_column=c('followT'),
    simplify=TRUE, advanced_dist=NULL)
}
\arguments{
  \item{group}{a character vector of the names of each group (e.g., \code{c('treatment','control')}).}
  \item{strata}{a character vector of the names of strata in groups (e.g., \code{c('young','old')}).}
  \item{allocation}{the relative ratio of sample size in each subgroup (\code{group*strata}). See details. The value will be recycled if the length is less than needed.}
  \item{event_lambda}{the hazard rate of the primary endpoint (event). See details. The value will be recycled if the length is less than needed.}
  \item{drop_rate}{(optional) the drop-out rate (patients/month). Not hazard rate. See details. The value will be recycled if the length is less than needed.}
  \item{death_lambda}{(optional) the hazard rate of death. The value will be recycled if the length is less than needed.}
  \item{n_rand}{(required when \code{rand_rate=NULL}) a vector of the number of randomization each month.}
  \item{rand_rate}{(required when \code{n_rand=NULL}) the randomization rate (patients/month).}
  \item{total_sample}{(required when \code{n_rand=NULL}) total scheduled sample size.}
  \item{add_column}{request additional columns of the returned data frame. \cr
Valid options are:
    \itemize{
      \item \code{'eventT_abs'}: absolute event time from the beginning of the trial (=eventT+randT)
      \item \code{'dropT_abs'}: absolute drop-out time from the beginning of the trial (=dropT+randT)
      \item \code{'deathT_abs'}: absolute death time from the beginning of the trial (=deathT+randT)
      \item \code{'censor'}: censoring (drop-out or death) indicator
      \item \code{'event'}: event indicator
      \item \code{'censor_reason'}: censoring reason ('drop_out','death','never_event'(eventT=inf))
      \item \code{'followT'}: follow-up time (true observed time) from randT
      \item \code{'followT_abs'}: absolute follow-up time from the beginning of the trial (=followT+randT)
    }}
  \item{simplify}{whether drop unused columns (e.g., the group variable when there is only one group). See details.}
  \item{advanced_dist}{use user-specified distributions for event, drop-out and death. A list containing random generation functions. See details and examples. }
}

\value{
  A data frame containing the some of these columns:
  \item{ID}{subject ID}
  \item{group}{group indicator}
  \item{strata}{stratum indicator}
  \item{randT}{randomization time (from the beginning of the trial)}
  \item{eventT}{event time (from \code{randT})}
  \item{eventT_abs}{event time (from the beginning of the trial)}
  \item{dropT}{drop-out time (from \code{randT})}
  \item{dropT_abs}{drop-out time (from the beginning of the trial)}
  \item{deathT}{death time (from \code{randT})}
  \item{deathT_abs}{death time (from the beginning of the trial)}
  \item{censor}{censoring (drop-out or death) indicator}
  \item{censor_reason}{censoring reason ('drop_out','death','never_event'(followT=inf))}
  \item{event}{event indicator}
  \item{followT}{follow-up time / observed time (from \code{randT})}
  \item{followT_abs}{follow-up time / observed time (from the beginning of the trial)}
}
\details{
See webpage \url{https://zjph602xtc.github.io/PWEXP/} for a diagram illustration of the relationship between returned variables.

The total number of subgroups will be '# treatment groups' * '# strata'. The \code{strata} variable will be distributed into each treatment group. For example, if \code{group = c('trt','placebo')}, \code{strata=c('A','B','C')}, then there will be 6 subgroups: trt+A, trt+B, trt+C, placebo+A, placebo+B, placebo+C. The lengths of \code{allocation}, \code{event_lambda}, \code{drop_rate}, \code{death_lambda} should be 6 as well. Note that the values will be recycled for these variables. For example, if \code{allocation=c(1,2,3)}, then the proportion of 6 subgroups is actually 1:2:3:1:2:3, which means 1:1 ratio for groups, 1:2:3 ratio in each stratum.

The \code{event_lambda} (\eqn{\lambda}) is the hazard rate of the interested events. The density function of events is \eqn{f(t)=\lambda e^{-\lambda*t}}. Similarly, the \code{death_lambda} is the hazard rate of death.

The \code{drop_rate} is the probability of drop-out at \eqn{t=1}, which means the hazard rate of drop-out is \eqn{-log(1-drop_rate)} (or say, \code{drop_rate}=\eqn{1-e^{-hazard rate}}.

When \code{simplify=TRUE}, these columns will NOT be included:
    \itemize{
      \item \code{group} when only one group is specified
      \item \code{strata} when only one stratum is specified
      \item \code{eventT} when \code{event_lambda=NA}
      \item \code{dropT} when \code{drop_rate=NA}
      \item \code{deathT}  when \code{death_lambda=NA}
    }

\code{advanced_dist} is used to define non-exponential distributions for event, drop-out or death. It is a list containing at least one of the elements: \code{event_dist}, \code{drop_dist}, \code{death_dist}. Each element has random generation functions for each subgroups. For example, \code{advanced_dist=list(event_dist=c(function1, function2), drop_dist=c(function3, function4))}. Here \code{function1}, \code{function3} are the event, drop-out generation function for the first subgroup; \code{function2}, \code{function4} for the second. If there is a third subgroup, \code{function1}, \code{function3} will be reused.
Each data generation function (\code{functionX}) is a function with only one input argument \code{n} (sample size). If any of the \code{event_dist}, \code{drop_dist}, \code{death_dist} is missing, then we search for \code{event_lambda}, \code{drop_rate}, \code{death_lambda} to generate a exp distribution; if they are also missing, then corresponding variable will not be generated
.}


\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{
  \code{\link{rpwexp}},
  \code{\link{rpwexp_conditional}}
}

\note{
\code{event_lambda}, \code{drop_rate}, \code{death_lambda} can be 0, which means the corresponding subgroup will have an Inf value for each variable.
}
\examples{
# Two groups with two strata. In the treatment group, there is a treatment
# sensitive stratum and a non-sensitive stratum. In the placebo group, all
# subjects are the same. Treatment:place=1:2. Drop rate=1\% only in treatment group.
dat <- simdata(group=c('trt', 'place'), strata = c('sensitive','non-sensitive'),
               allocation = c(1,1,2,2), rand_rate = 20, total_sample = 1000,
               event_lambda = c(0.1, 0.2, 0.01, 0.01),
               drop_rate = c(0.01, 0.01, 0, 0))
# randomized subjects
table(dat$group,dat$strata)
# randomization curve
plot(sort(dat$randT), 1:1000, xlab='time', ylab='randomized subjects')
# event time in treatment group
plot(ecdf(dat$eventT[dat$group=='trt' & dat$strata=='sensitive']))
lines(ecdf(dat$eventT[dat$group=='trt' & dat$strata=='non-sensitive']), col='red')


# One group. Event follows a piecewise exponential distribution; drop-out follows
# a Weibull; death follows a exponential.
dist_trt <- function(n)rpwexp(n, rate=c(0.01, 0.05, 0.01), breakpoint = c(30,60))
dist_placebo <- function(n)rpwexp(n, rate=c(0.01, 0.005), breakpoint = c(50))
dat <- simdata(group = c('trt','placebo'), n_rand = c(rep(10,50),rep(20,10)),
               death_lambda = 0.01,
               advanced_dist = list(event_dist=c(dist_trt, dist_placebo),
                                    drop_dist=function(n)rweibull(n,3,40)))
# randomized subjects
table(dat$group)
# randomization curve
plot(sort(dat$randT), 1:700, xlab='time', ylab='randomized subjects')
# event time in both groups
plot(ecdf(dat$eventT[dat$group=='trt']), xlim=c(0,100))
lines(ecdf(dat$eventT[dat$group=='placebo']), col='red')
# drop-out time
plot(ecdf(dat$dropT), xlim=c(0,100))


# mixture cure distribution, 20\% of the subject are cured and will not have events
dat <- simdata(strata=c('cure','non-cure'), allocation=c(20,80),
        event_lambda=c(0, 0.38), n_rand = rep(20,30),
        add_column = c('eventT_abs', 'censor', 'event',
                       'censor_reason', 'followT', 'followT_abs'))
}










