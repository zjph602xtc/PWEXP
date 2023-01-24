\name{plot_event}
\title{Plot Cumulative Event Curve}
\description{
  Plot cumulative event curve with right censoring data.
}
\usage{
plot_event(time, abs_time=TRUE, event, additional_event=0,
         add=FALSE, plot=TRUE, xyswitch=FALSE, ...)
}
\arguments{
 \item{time}{observed/follow-up time from the first day of the trial or individual time from enrollment.}
 \item{abs_time}{logical; if TRUE, \code{time} is the time from the first day of the trial. }
 \item{event}{the status indicator, 0=censor, 1=event. }
 \item{additional_event}{adding the cumulative number of events by a constant number from the beginning. }
 \item{add}{logical; if TRUE add lines to current plot.}
 \item{plot}{logical; if FALSE, do not plot any lines, but return the line data}
 \item{xyswitch}{logical; if TRUE, x-axis will be cumulative number of events and y will be the time.}
 \item{...}{pass over to \code{\link[base]{plot}.}}
}
\details{
A convenient function to calculate and plot the cumulative number of events.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\value{
 A data frame containing these columns:
  \item{time}{sorted time}
  \item{n_event}{cumulative number of events}
}
\seealso{\code{\link{plot_km}}}
