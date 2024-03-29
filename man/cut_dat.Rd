\name{cut_dat}
\alias{cut_dat}
\title{Cut Data before a Specified Time}
\description{
  Take a subset of a dataset by constraining the randomization time <= cut time. Additionally, it updates the follow-up time, censor/event indicator, censor reason, accordingly.
}
\usage{
cut_dat(cut, data, var_randT=NULL, var_followT=NULL, var_followT_abs=NULL,
      var_censor=NULL, var_event=NULL, var_censor_reason='status_at_end')
}
\arguments{
 \item{cut}{cut time (from the beginning of the trial); only rows with randomization time <= \code{cut} will be kept.}
 \item{data}{a data frame.}
 \item{var_randT}{character; the variable name of randomization time. If missing, then the randomization time will be treated as 0 and NO subjects will be filtered by \code{cut} time.}
 \item{var_followT}{character; the variable name of follow-up time (from randomization)}
 \item{var_followT_abs}{character; the variable name of follow-up time (from the beginning of the trial)}
 \item{var_censor}{character; the variable name of censoring (drop-out or death) indicator (1=censor, 0=event)}
 \item{var_event}{character; the variable name of event indicator (1=event, 0=censor)}
 \item{var_censor_reason}{character; the variable name of censoring reason (character variable). This variable will be created, if \code{data} does not contain it.}
}
\details{
We first filter rows that randomization time is equal to or less then \code{cut} time. Then we modify these columns (if provided):
\itemize{
      \item \code{var_followT:} change values to (\code{cut} - randomization time) if (follow-up time + randomization time) > \code{cut}
      \item \code{var_followT_abs:} change values to \code{cut} if (follow-up time from beginning) > \code{cut}
      \item \code{var_censor:} change values to 1 if (follow-up time from beginning) > \code{cut}
      \item \code{var_event:} change values to 0 if (follow-up time from beginning) > \code{cut}
      \item \code{var_censor_reason:} change values to 'cut' if (follow-up time from beginning) > \code{cut}
    }
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\value{
 A subset data frame with the same columns as \code{data}.

 \code{var_censor_reason} is the only variable that is allowed to be absent in \code{data}. The function will create this variable in the returned data frame and set values 'cut' to the subjects whose (follow-up time from beginning) > \code{cut}.
}


\examples{
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(rand_rate = 20, total_sample = 1000, drop_rate = 0.03,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))
cut <- quantile(dat$randT, 0.8)
train <- cut_dat(var_randT = 'randT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')
}


\note{
The original dataset \code{data} will NOT be modified.
}
