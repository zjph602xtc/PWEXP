\name{predict}
\alias{predict.pwexp.fit}
\alias{predict.boot.pwexp.fit}
\title{Predict Events for Piecewise Exponential Model}
\description{
 Obtians event prediction and (optionally) confidence interval from a piecewise exponential model.
}
\usage{
\method{predict}{pwexp.fit}(event_model, cut_indicator, cut_time, censor_model,
      n_each=100, future_enroll, seed=1818, plot=F)
\method{predict}{boot.pwexp.fit}(event_model_boot, cut_indicator, cut_time,
      censor_model_boot, n_each=10, future_enroll, seed=1818, plot=F)
}
\arguments{
  \item{event_model}{an object of class \code{pwexp.fit} returned by the \code{\link{pwexp.fit}} function. It is the event model for the primary endpoint. }
    \item{event_model_boot}{an object of class \code{boot.pwexp.fit} returned by the \code{\link{boot.pwexp.fit}} function. It is the event model with bootstrapping for the primary endpoint. }
  \item{cut_indicator}{(optional) A vector indicates which subject is censored due to the end of the trial. The length of the vector is the number of rows of the data used in \code{event_model}/\code{event_model_boot}. Value 0 means the subject had event or drop-out or death before the end of the trial; 1 means the subject didn't have any of these. See details. }
  \item{cut_time}{the weight of each subject. It should be in the same length of \code{train.y}.}
  \item{censor_model_boot}{if \code{FALSE}, any subjects in the training data and the validation data (if exist) with zero (or negative) weights will be removed.}
  \item{n_each}{if \code{TRUE}, we prefit the model with partitioned training data using the first set of parameters in \code{range}. If fails (i.e., too many zero weight subjects in the partitioned training data), we re-partition the data and re-try the model for up to 10 times. This is useful when \code{use_zero_weight=TRUE} and there many zero weights subjects in the data. }
  \item{future_enroll}{data, if a formula interface is used. Ignored, if predictor matrix and response are supplied directly.}
  \item{seed}{a random seed.}
  \item{plot}{}
}
\value{
  For \code{tune_wsvm}, an object of class \code{tune_wsvm}, including the components:
  \item{best.parameters}{a 1 x k data frame, k number of parameters.}
  \item{best.performance}{best achieved performance.}
  \item{performances}{if requested, a data frame of all parameter
    combinations along with the corresponding performance results.}
  \item{train.ind}{list of index vectors used for splits into
    training and validation sets.}
  \item{best.model}{if requested, the model trained on the complete training data
    using the best parameter combination.}

  \code{best.tune_wsvm()} returns the best model detected by \code{tune_wsvm}.
}
\details{
\code{cut_indicator} indicates the status of each subject in the \code{event_model}/\code{event_model_boot} model at the end of the trial. Value 1 means the subject didn't have events, drop-out or death at the end of the trial (or say, the subject was censored due to the end of the trial). When \code{cut_indicator} is NOT provided, we assign value 1 to the subject who didn't have event (or drop-out, or death) in both  \code{event_model}/\code{event_model_boot} and \code{censor_model}/\code{censor_model_boot} models.
}
\author{
  David Meyer \cr
  Modified by Tianchen Xu \email{tx2155@columbia.edu}
}
\seealso{\code{\link{tune.control}}, \code{\link{plot.tune_wsvm}}}
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
                     cut_time = cut,  censor_model_boot=fit_res_censor_boot,
                     future_enroll=list(enroll_rate=20, total_sample=NROW(dat)-NROW(train)))

}
