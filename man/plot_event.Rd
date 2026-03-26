\name{plot_event}
\alias{plot_event}
\alias{plot_event.default}
\alias{plot_event.predict.pwexpm}
\alias{plot_event.predict.boot.pwexpm}
\alias{plot.predict.pwexpm}
\alias{plot.predict.boot.pwexpm}
\title{Plot Cumulative Event Curve with Survival Data or Predicted PWE Model}
\description{
  Plot cumulative event curve with right censoring data.
}
\usage{
\method{plot_event}{default}(time, event, abs_time=TRUE, additional_event=0,
         add=FALSE, plot=TRUE, xyswitch=FALSE, ...)
\method{plot_event}{predict.pwexpm}(time, abs_time=TRUE, add=TRUE, plot=TRUE,
         xyswitch=FALSE, eval_at=NULL, ...)
\method{plot}{predict.pwexpm}(x, ...)
\method{plot_event}{predict.boot.pwexpm}(time, abs_time=TRUE, alpha=0.1, type='confidence',
         add=TRUE, plot=TRUE, xyswitch=FALSE, eval_at=NULL,
         show_CI=TRUE, CI_par=NULL, ...)
\method{plot}{predict.boot.pwexpm}(x, ...)
}
\arguments{
 \item{time}{observed/follow-up time from individual randomization time (\code{abs_time=FALSE}) or from the first subject randomization time (\code{abs_time=TRUE}); or a predicted object from \code{\link{predict.pwexpm}}, or a predicted object with bootstrapping  from \code{\link{predict.boot.pwexpm}}.}
 \item{x}{the same as \code{time}.}
 \item{abs_time}{logical; if TRUE, \code{time} is the time from first randomization of the trial. if FALSE, \code{time} is the time from the randomization of each subject. }
 \item{event}{the status indicator, 0=censor, 1=event. Other choices are TRUE/FALSE (TRUE = event).}
 \item{additional_event}{adding the cumulative number of events by a constant number from the beginning. }
 \item{add}{logical; if TRUE add lines to current plot.}
 \item{plot}{logical; if FALSE, do not plot any lines, but return the line data}
 \item{xyswitch}{logical; if TRUE, x-axis will be cumulative number of events and y will be the time.}
 \item{eval_at}{a vector of the time (when \code{xyswitch=FALSE}) or the number of events (when \code{xyswitch=TRUE}) that you want to make prediction on. }
 \item{alpha}{the significance level of the confidence interval.}
 \item{type}{the type of prediction required. The default \code{confidence} returns the confidence interval without random error; the alternative \code{predictive} returns the predictive interval.}
 \item{show_CI}{logical; if TRUE add confidence interval of the estimated event curve. }
 \item{CI_par}{a list of parameters to control the apperance of lines of confidence intervals. The values pass to \code{\link[graphics]{lines}}.}
 \item{...}{other arguments (e.g., \code{lwd}, etc.) are passed over to \code{\link[graphics]{plot.default}.}}
}
\details{
A convenient function to calculate and plot the cumulative number of events.

Parameters in \code{...} are passed to \code{\link[graphics]{plot.default}} function to control the appearance of the event curve; parameters in \code{CI_par} are passed to \code{\link[graphics]{lines}} function to control the appearance of confidence intervals. See examples for usage.

By default, \code{plot_event} plots a data frame in a new figure; and plots a predicted model in existing figure.

\code{plot} method for classes \code{predict.pwexpm} and \code{predict.boot.pwexpm} behaves the same as \code{plot_event}, except that it plots a predicted model in a new figure by default.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\value{
When \code{xyswith=FALSE}, the function returns a data frame containing these columns:
  \item{time}{sorted time at \code{eval_at}}
  \item{n_event}{predicted cumulative number of events}
  \item{\code{(alpha/2)} n_event}{(for \code{predict.boot.pwexpm}) lower alpha level predicted cumulative number of events}
  \item{\code{1-(alpha/2)} n_event}{(for \code{predict.boot.pwexpm}) upper alpha level predicted cumulative number of events}
When \code{xyswith=TRUE}, the function returns a data frame containing these columns:
  \item{n_event}{sorted cumulative number of events at \code{eval_at}}
  \item{time}{predicted required time}
  \item{\code{(alpha/2)} time}{(for \code{predict.boot.pwexpm}) lower alpha level predicted required time}
  \item{\code{1-(alpha/2)} time}{(for \code{predict.boot.pwexpm}) upper alpha level predicted required time}
}

\examples{
set.seed(1818)
event_dist <- function(n)rpwexpm(n, rate = c(0.1, 0.2), breakpoint = 14)
dat <- simdata(rand_rate = 20, drop_rate = 0.03, total_sample = 500,
               advanced_dist = list(event_dist=event_dist),
               add_column = c('censor_reason','event','followT','followT_abs'))
cut <- quantile(dat$randT, 0.8)
train <- cut_dat(var_randT = 'randT', cut = cut, data = dat,
                 var_followT = 'followT', var_followT_abs = 'followT_abs',
                 var_event = 'event', var_censor_reason = 'censor_reason')

fit_res3 <- pwexpm(Surv(followT, event), data=train, nbreak = 1)
fit_res_boot <- boot.pwexpm(fit_res3, nsim = 8) # here nsim=8 is for demo purpose,
# pls increase it in practice

drop_indicator <- ifelse(train$censor_reason=='drop_out' & !is.na(train$censor_reason),1,0)
fit_res_censor <- pwexpm_fit(train$followT, drop_indicator, nbreak = 0)
fit_res_censor_boot <- boot.pwexpm(fit_res_censor, nsim = 8)

cut_indicator <- train$censor_reason=='cut'
cut_indicator[is.na(cut_indicator)] <- 0

predicted_boot <- predict(fit_res_boot, cut_indicator = cut_indicator,
                          analysis_time = cut,  censor_model=fit_res_censor_boot,
                          future_rand=list(rand_rate=20, total_sample=NROW(dat)-NROW(train)))

plot_event(train$followT_abs, train$event, xlim=c(0,69), ylim=c(0,500))
plot_event(predicted_boot, eval_at = seq(40,90,5), CI_par = list(lty=3, lwd=2))

plot_event(train$followT_abs, train$event, xyswitch = TRUE, ylim=c(0,69), xlim=c(0,400))
plot_event(predicted_boot, xyswitch = TRUE,  eval_at = seq(250,400,50))
}

\seealso{\code{\link{plot_survival}}}
