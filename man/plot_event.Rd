\name{plot_event}
\alias{plot_event}
\alias{plot_event.default}
\alias{plot_event.predict.pwexp.fit}
\alias{plot_event.predict.boot.pwexp.fit}
\title{Plot Cumulative Event Curve}
\description{
  Plot cumulative event curve with right censoring data.
}
\usage{
\method{plot_event}{default}(time, abs_time=TRUE, event, additional_event=0,
         add=FALSE, plot=TRUE, xyswitch=FALSE, ...)
\method{plot_event}{predict.pwexp.fit}()
\method{plot_event}{predict.boot.pwexp.fit}()
}
\arguments{
 \item{time}{observed/follow-up time from  the enrollment time or from the beginning of the trial.}
 \item{abs_time}{logical; if TRUE, \code{time} is the time from the beginning of the trial. }
 \item{event}{the status indicator, 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event).}
 \item{additional_event}{adding the cumulative number of events by a constant number from the beginning. }
 \item{add}{logical; if TRUE add lines to current plot.}
 \item{plot}{logical; if FALSE, do not plot any lines, but return the line data}
 \item{xyswitch}{logical; if TRUE, x-axis will be cumulative number of events and y will be the time.}
 \item{...}{other arguments (e.g., \code{lwd}, etc.) are passed over to \code{\link[base]{plot}.}}
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

\examples{
# plot number of events for dataset
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(enroll_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))

plot_event(dat$followT_abs, dat$event, abs_time=T, xlim=c(0,23), ylim=c(0,290))
plot_event(dat$followT, dat$event, abs_time=F, xlim=c(0,23), ylim=c(0,1000))

# plot predicted number of events
cut <- quantile(dat$enrollT, 0.8)
train <- cut_dat(var_enrollT = 'enrollT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')
fit_res3 <- pwexp.fit(train$followT, train$event, nbreak = 2)
drop_indicator <-  ifelse(train$censor_reason=='drop_out',1,0)
drop_indicator[is.na(drop_indicator)] <- 0
fit_res_censor <- pwexp.fit(train$followT, drop_indicator, nbreak = 0)
cut_indicator <- train$censor_reason=='cut'
cut_indicator[is.na(cut_indicator)] <- 0
predicted <- predict(event_model=fit_res3, cut_indicator = cut_indicator,
                     cut_time = cut,  censor_model=fit_res_censor,
                     future_enroll=list(enroll_rate=20, total_sample=NROW(dat)-NROW(train)))
plot_event(predicted, eval_at = 40:90, plot = T, add = F, xlim=c(40,90)
}


\seealso{\code{\link{plot_survival}}}
