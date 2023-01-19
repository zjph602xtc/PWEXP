\name{plot_km}
\title{Plot KM Curve with Right Censoring}
\description{
  Plot KM Curve.
}
\usage{
plot_km(time, event, conf.int=FALSE, mark.time=TRUE, lwd=2,
      xlab='Follow-up time', ylab='Survival function', ...)
}
\arguments{
 \item{time}{observed time from enrollment. For right censored data, this is the follow-up time. }
 \item{event}{the status indicator, normally 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event) or 1/2 (2=event). }
 \item{...}{pass over to \code{\link[survival]{plot.survfit}.}}
}
\details{
A wrapper of \code{\link[survival]{plot.survfit}} function to plot right censoring data.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}

\seealso{\code{\link{plot_event}}
