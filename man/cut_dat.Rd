\name{cut_dat}
\title{Cut Data by Enrollment Time}
\description{
  Take a subset of a dataset by constraining the enrollment time <= cut time. Additionally, it updates the follow-up time, censor/event indicator, censor reason, accordingly.
}
\usage{
cut_dat(var_enrollT, cut, data, var_followT=NULL, var_followT_abs=NULL,
      var_censor=NULL, var_event=NULL, var_censor_reason=NULL)
}
\arguments{
 \item{var_enrollT}{character; the variable name of enrollment time}
 \item{cut}{cut time (from the beginning of the trial); only rows with enrollment time <= \code{cut} will be kept.}
 \item{data}{a data frame.}
 \item{var_followT}{character; the variable name of follow-up time (from enrollment)}
 \item{var_followT_abs}{character; the variable name of follow-up time (from the beginning of the trial)}
 \item{var_censor}{character; the variable name of censoring (drop-out or death) indicator (1=censor, 0=event)}
 \item{var_event}{character; the variable name of event indicator (1=event, 0=censor)}
 \item{var_censor_reason}{character; the variable name of censoring reason (character variable)}
}
\details{
We first filter rows that enrollment time is equal to or less then \code{cut} time. Then we modify these columns (if provided):
\itemize{
      \item \code{var_followT:} change values to (\code{cut} - enrollment time) if (follow-up time + enrollment time) > \code{cut}
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
}

\note{
The original dataset \code{data} will NOT be modified.
}

\seealso{\code{\link{tran_censor}}}
