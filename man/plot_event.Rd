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
\method{plot_event}{predict.pwexp.fit}(predict_model, abs_time=TRUE, add=TRUE, plot=TRUE,
         xyswitch=FALSE, eval_at=NULL, ...)
\method{plot_event}{predict.boot.pwexp.fit}(predict_model_boot, abs_time=TRUE, alpha=0.05, add=TRUE,
         plot=TRUE, xyswitch=FALSE, eval_at=NULL, show_CI=TRUE, CI_par=NULL, ...)
}
\arguments{
 \item{time}{observed/follow-up time from  the enrollment time or from the beginning of the trial.}
 \item{abs_time}{logical; if TRUE, \code{time} is the time from the beginning of the trial. }
 \item{event}{the status indicator, 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event).}
 \item{predict_model}{a predicted object of class inheriting from \code{\link{predict.pwexp.fit}}.}
 \item{predict_model_boot}{a predicted object with bootstrapping of class inheriting from \code{\link{predict.boot.pwexp.fit}}.}
 \item{additional_event}{adding the cumulative number of events by a constant number from the beginning. }
 \item{add}{logical; if TRUE add lines to current plot.}
 \item{plot}{logical; if FALSE, do not plot any lines, but return the line data}
 \item{xyswitch}{logical; if TRUE, x-axis will be cumulative number of events and y will be the time.}
 \item{eval_at}{a vector of the time (when \code{xyswitch=FALSE}) or the number of events (when \code{xyswitch=TRUE}) that you want to make prediction on. }
 \item{alpha}{the significance level of the confidence interval.}
 \item{show_CI}{logical; if TRUE add confidence interval of the estimated event curve. }
 \item{CI_par}{a list of parameters to control the apperance of lines of confidence intervals. The values pass to \code{\link[graphics]{lines}}.}
 \item{...}{other arguments (e.g., \code{lwd}, etc.) are passed over to \code{\link[base]{plot}.}}
}
\details{
A convenient function to calculate and plot the cumulative number of events.

Parameters in \code{...} are passed to \code{\link[base]{plot}} function to control the appearance of the event curve; parameters in \code{CI_par} are passed to \code{\link[graphics]{lines}} function to control the appearance of confidence intervals. See examples for usage.

By default, \code{plot_event} plots a data frame in a new figure; and plots a predicted model in existing figure.
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
set.seed(1818)
event_dist <- function(n)rpwexp(n, rate = c(0.1, 0.01, 0.2), breakpoint =  c(5,14))
dat <- simdata(enroll_rate = 20, drop_rate = 0.03,  total_sample = 1000,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))
cut <- quantile(dat$enrollT, 0.8)
train <- cut_dat(var_enrollT = 'enrollT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')

fit_res3 <- pwexp.fit(train$followT, train$event, nbreak = 2)
fit_res_boot <- boot.pwexp.fit(fit_res3, nsim = 100)

drop_indicator <-  ifelse(train$censor_reason=='drop_out' & !is.na(train$censor_reason),1,0)
fit_res_censor <- pwexp.fit(train$followT, drop_indicator, nbreak = 0)
fit_res_censor_boot <- boot.pwexp.fit(fit_res_censor)

cut_indicator <- train$censor_reason=='cut'
cut_indicator[is.na(cut_indicator)] <- 0

predicted_boot <- predict(event_model_boot=fit_res_boot, cut_indicator = cut_indicator,
                     analysis_time = cut,  censor_model_boot=fit_res_censor_boot,
                     future_enroll=list(enroll_rate=20, total_sample=NROW(dat)-NROW(train)))

plot_event(dat$followT_abs, abs_time = T, event=dat$event, ylim=c(0,800))
plot_event(predicted_boot, eval_at = seq(40,90,5), CI_par = list(lty=3, lwd=2))

plot_event(dat$followT_abs, xlim=c(0,800), event=dat$event, xyswitch = T)
plot_event(predicted_boot, xyswitch = T, eval_at = seq(600,900,50))
}


\seealso{\code{\link{plot_survival}}}
