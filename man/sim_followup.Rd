\name{sim_followup}
\alias{sim_followup}
\alias{print.sim_followup}
\title{Estimate Follow-up Time and Number of Events by Simulation}
\description{
\code{sim_follwup} is used to estimate follow-up time and number of events (given calander time, or number of randomized samples, or number of events).
}
\usage{
sim_followup(at, type = 'calander', group="Group 1", strata='Strata 1',
            allocation=1, event_lambda=NA, drop_rate=NA, death_lambda=NA,
            n_rand=NULL, rand_rate=NULL, total_sample=NULL, extra_follow=0,
            by_group=FALSE, by_strata=FALSE, advanced_dist=NULL,
            stat=c(mean, median, sum), follow_up_endpoint=c('death', 'drop_out',
            'cut'), count_in_extra_follow=FALSE, count_insufficient_event=FALSE,
            start_date=NULL, rep=300, seed=1818)
}
\arguments{
  \item{at}{specify a vector of occasions. When \code{type='calander'}, \code{at} is the time from fisrt randomization; when \code{type='event'}, \code{at} is the number of accumulated events; when \code{type='sample'}, \code{at} is  the number of randomized samples. }
  \item{type}{specify the type of \code{at}. Must be \code{'calander'}, \code{event} or \code{sample}.}
  \item{group}{a character vector of the names of each group (e.g., \code{c('treatment','control')}). See \code{\link{simdata}}.}
  \item{strata}{a character vector of the names of strata in groups (e.g., \code{c('young','old')}). See \code{\link{simdata}}.}
  \item{allocation}{the relative ratio of sample size in each subgroup (\code{group*strata}).  The value will be recycled if the length is less than needed. See \code{\link{simdata}}.}
  \item{event_lambda}{the hazard rate of the primary endpoint (event). The value will be recycled if the length is less than needed. See \code{\link{simdata}}.}
  \item{drop_rate}{(optional) the drop-out rate (patients/month). Not hazard rate. The value will be recycled if the length is less than needed. See \code{\link{simdata}}.}
  \item{death_lambda}{(optional) the hazard rate of death. The value will be recycled if the length is less than needed. See \code{\link{simdata}}.}
  \item{n_rand}{(required when \code{rand_rate=NULL}) a vector of the number of randomization each month; can be non-integers. See \code{\link{simdata}}.}
  \item{rand_rate}{(required when \code{n_rand=NULL}) the randomization rate (patients/month; can be non-integer). See \code{\link{simdata}}.}
  \item{total_sample}{(required when \code{n_rand=NULL}) total scheduled sample size. See \code{\link{simdata}}.}
  \item{extra_follow}{delay the analysis time by extra time (\code{extra_follow}) after the time specified by \code{at}. See details. }
  \item{by_group}{logical; if TRUE, also return results by each group. }
  \item{by_strata}{logical; if TRUE, also return results by each stratum. }
  \item{advanced_dist}{use user-specified distributions for event, drop-out and death. A list containing random generation functions. See details and examples in \code{\link{simdata}}.}
  \item{stat}{a vector of functions to summarize the follow-up time. See example. }
  \item{follow_up_endpoint}{Which endpoints can be regarded as the end of follow-up. Choose from 'death', 'drop_out', 'cut' (censored at the end of the trial) or 'event'.'}
  \item{count_in_extra_follow}{logical; whether to count subjects who are randomized after the time spcified by \code{at} but before the time specified by \code{at} + \code{extra_follow}.}
  \item{count_insufficient_event}{logical; only affects the result when \code{type='event'}. If TRUE, for samples that cannot achieve required number of events, the last follow-up time is the analysis time. If FALSE, these samples will be dropped. }
  \item{start_date}{the start date of the first randomization; in the format: "2000-01-30"}
  \item{rep}{number simulated iterations. }
  \item{seed}{a random seed. Do not set seed if \code{seed=NULL}.}
}

\value{
  A list contains some of the data frames:
  \item{T_all}{Simulation result for the whole population}
  \item{T_by_group}{Simulation result by arms}
  \item{T_by_strata}{Simulation result by strata}
  \item{T_by_group_strata}{Simulation result by arms and strata}

  Each data frame contains some of these columns:
  \item{at}{When \code{type='calander'}, \code{at} is the time from fisrt randomization; when \code{type='event'}, \code{at} is the number of accumulated events; when \code{type='sample'}, \code{at} is  the number of randomized samples}
  \item{group}{group indicator}
  \item{strata}{stratum indicator}
  \item{analysis_time}{time from the first randomization}
  \item{analysis_time_c}{calendar time from the first randomization}
  \item{event}{number of events}
  \item{subjects}{number of enrolled subjects}
  \item{functions in \code{stat}}{function names specified in the \code{stat} argument}
}
\details{
See the help document of \code{\link{simdata}} for most arguments details.

When \code{type='calander'}, the function estimates the follow-up time and number of events at time \code{at} plus \code{extra_follow}; when \code{type='event'},  the function estimates these at the time when total number of events is \code{at} plus time \code{extra_follow}; when \code{type='sample'}, the function estimates these at the time when total number of randomized subjects is \code{at} plus time \code{extra_follow}.

The \code{stat} specifies a vector of user defined functions. Each of them must take a vector of individual follow-up time as input and return a single summary value. See example.
}


\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{simdata}}}

\note{
\code{event_lambda}, \code{drop_rate}, \code{death_lambda} can be 0, which means the corresponding subgroup will have no events, drop-out, or deaths.
}
\examples{
# Two groups. Treatment:place=1:2. Drop rate=3\%/month. Hazard ratio=0.7.

# define the piecewiese exponential event generation function
myevent_dist_trt <- function(n)rpwexpm(n, rate=c(0.1, 0.01, 0.2)*0.7, breakpoint=c(5,14))
myevent_dist_con <- function(n)rpwexpm(n, rate=c(0.1, 0.01, 0.2), breakpoint=c(5,14))

# user defined summary function, the proportion of subjects that follow more than 12 month
prop_12 <- function(x)mean(x >= 12)

# estimate the event curve or timeline:
# (here rep=60 is for demo purpose only, please increase this value in practice!)
event_curve <- sim_followup(at=seq(20,90,10), type = 'calendar', group = c('trt','con'),
             rand_rate = 20, total_sample = 1000, drop_rate = 0.03, allocation = 1:2,
             advanced_dist = list(event_dist=c(myevent_dist_trt, myevent_dist_con)),
             by_group = TRUE, stat = c(median, mean, prop_12), start_date = "2020-01-01",
             rep=60)
time_curve <- sim_followup(at=seq(200,600,100), type = 'event', group = c('trt','con'),
             rand_rate = 20, total_sample = 1000, drop_rate = 0.03, allocation = 1:2,
             advanced_dist = list(event_dist=c(myevent_dist_trt, myevent_dist_con)),
             stat = c(median, mean, prop_12), start_date = "2020-01-01", rep=60)
# plot event curve or timeline
plot(event_curve$T_all$analysis_time_c, event_curve$T_all$event, xlab='Time',
    ylab='Number of events', type='b')
plot(time_curve$T_all$event, time_curve$T_all$analysis_time_c, xlab='Number of
    events', ylab='Time', type='b')
}










