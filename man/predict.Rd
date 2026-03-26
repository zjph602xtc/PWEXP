\name{predict}
\alias{predict.pwexpm}
\alias{predict.boot.pwexpm}
\alias{summary.predict.pwexpm}
\alias{summary.predict.boot.pwexpm}
\title{Predict Events for Piecewise Exponential (PWE) Model}
\description{
Obtain event prediction and (optionally) confidence interval from a piecewise exponential model.
}
\usage{
\method{predict}{pwexpm}(object, cut_indicator=NULL, analysis_time, censor_model=NULL,
      n_each=100, future_rand=NULL, seed=1818, ...)
\method{predict}{boot.pwexpm}(object, cut_indicator=NULL, analysis_time, censor_model=NULL,
      n_each=10, future_rand=NULL, seed=1818, ...)
}
\arguments{
  \item{object}{a \code{\link{pwexpm}} or \code{\link{boot.pwexpm}} object. It is the event model for the primary endpoint. }
  \item{cut_indicator}{(optional) A vector indicates which subject is censored due to the end of the trial. The length of the vector is the number of rows of the data used in \code{event_model}/\code{event_model_boot}. Value 0 means the subject had event or drop-out or death before the end of the trial; 1 means the subject didn't have any of these. See details. }
  \item{analysis_time}{the analysis time. This is the time length from the start of the trial to the time collecting data for the model. }
  \item{censor_model}{an object of class \code{\link{pwexpm}} (or \code{\link{boot.pwexpm}}). It is the censoring model for drop-out and death. }
  \item{n_each}{the number of iterations for each bootstrapping sample to obtain predicitive CI. Typically, a value of 10 to 100 should be enough. }
  \item{future_rand}{the randomization curve in the following times. Can be \code{NULL} if all subjects have been randomized. You can specify \bold{future} rand rate and \bold{future} total number of samples to be randomized by \code{list(rand_rate= ,total_sample= )} or specify the \bold{future} number of randomization each month by \code{list(n_rand= )}. See details.}
  \item{seed}{a random seed. Do not set seed if \code{seed=NULL}.}
  \item{...}{internal function reserved. }
}
\value{
A object of class "\code{predict.pwexpm}" or "\code{predict.boot.pwexpm}" is a list containing the following components:
  \item{event_fun}{ number of events vs. time curve function in each bootstrap.}
  \item{event_model}{the event model for the primary endpoint.}
  \item{censor_model}{the censoring model for drop-out and death.}
  \item{nsim}{the number of repeated bootstrapping. \code{nsim=1} for non-bootstrapped model.}
  \item{bootstrap}{a logical value indicating if the \code{event_model} is a bootstrapped model.}
  \item{para}{the parameters used to conduct the prediction procedure.}

  This returned object should be used in \code{\link{plot_event}} function for summarizing its result.
}
\details{
The prediction will have a confidence interval only if the event model and censor model are bootstrap models.

\code{cut_indicator} indicates the status of each subject in the \code{event_model}/\code{event_model_boot} model at the end of the trial. Value 1 means the subject didn't have events, drop-out or death at the end of the trial (or say, the subject was censored due to the end of the trial). When \code{cut_indicator} is NOT provided, we assign value 1 to the subject who didn't have event (or drop-out, or death) in both  \code{event_model}/\code{event_model_boot} and \code{censor_model}/\code{censor_model_boot} models.

\code{future_rand} is a list determining the parameter of randomization curve in the following times. For example, you specify randomization rate=10pt/month and total sample size=1000 by \code{list(rand_rate=10, total_sample=1000)} or specify the number of randomization each month (e.g., 10,15,30,30 in four months) by \code{list(n_rand=c(10,15,30,30))}.
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{\code{\link{plot_event}}}
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
