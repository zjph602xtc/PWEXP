\name{simdata}
\title{Simuate Survival Data}
\description{
\code{simdata} is used to simuate a clinical trial data with time-to-event endpoints.
}
\usage{
\method{wsvm}{default}(group="Group 1", strata='Strata 1', allocation=1,
    event_lambda=NA, drop_rate=NA, death_lambda=NA, n_enroll=NULL,
    enroll_rate=NULL, total_sample=NULL, add_column=NULL, simplify=TRUE,
    advanced_dist=NULL)
}
\arguments{
  \item{group}{group}
  \item{strata}{strata}
  \item{allocation}{the proportion of each subgroup. For example, when group=c('trt','placebo'), strata=c('A','B','C'), the length of allocation should be 6. Note the value will be recycled if the length is less than needed.}
  \item{event_lambda}{the hazard rate of event. The value will be recycled if the length is less than needed.}
  \item{drop_rate}{the drop-out rate (patient/month). The value will be recycled if the length is less than needed.}
  \item{death_lambda}{the hazard rate of death. The value will be recycled if the length is less than needed.}
  \item{n_enroll}{(required when enroll_rate=NULL) a vector contains the number of enrollment each month}
  \item{enroll_rate}{(required when n_enroll=NULL) the enrollment rate (patient/month)}
  \item{total_sample}{(required when enroll_rate is in use) total scheduled sample size}
  \item{add_column}{request additional columns \cr
Valid options are:
    \itemize{
      \item \code{eventT_abs}absolute event time from the beginning of the trial (=eventT+enrollT)
      \item \code{deathT_abs}absolute death time from the beginning of the trial (=deathT+enroolT)
      \item \code{censor} whether censored
      \item \code{event}whether having event
      \item \code{censor_reason}why censored ('censored','death','never_event'(eventT=inf))
      \item \code{followT}follow time (true observed time) from enrollT
      \item \code{followT_abs}absolute follow time from the beginning of the trial (=followT+enrollT)
    }}
  \item{simplify}{whether drop unused columns (i.e., the group variable when there is only one group)}
  \item{advanced_dist}{use other distributions instead of exponential. The value will be recycled if the length is less than needed. e.g., list(event_dist=c(function(n)rexp(n,4), function(n)rweibull(n,3,2), drop_dist=..., death_dist=...) Each distribution is a function with only one input argument n (sample size) If any of event_dist, drop_dist, death_dist is missing, then search for event_lambda, drop_rate, death_lambda; if also missing, then this variable will not be generated}
}

\value{
  An object of class \code{"wsvm"} containing the fitted model, including:
  \item{SV}{The resulting support vectors (possibly scaled).}
  \item{index}{The index of the resulting support vectors in the data
    matrix. Note that this index refers to the preprocessed data (after
    the possible effect of \code{na.omit} and \code{subset})}
  \item{coefs}{The corresponding coefficients times the training labels.}
  \item{rho}{The negative intercept.}
  \item{sigma}{In case of a probabilistic regression model, the scale
    parameter of the hypothesized (zero-mean) laplace distribution estimated by
    maximum likelihood.}
  \item{probA, probB}{numeric vectors of length k(k-1)/2, k number of
    classes, containing the parameters of the logistic distributions fitted to
    the decision values of the binary classifiers (1 / (1 + exp(a x + b))).}
}
\details{
total number of subgroups will be '# treatment groups'*'# strata'  strata variable will be distributed into each treatment group. For example,
group <- c('trt','placebo'), strata <- c('A','B','C'), then there will be six subgroups:
 trt+A, trt+B, trt+C, placebo+A, placebo+B, placebo+C
}
\note{
Data are scaled internally, usually yielding better results.

Parameters of SVM-models usually \emph{must} be tuned to yield sensible results!
}
\references{
  \itemize{
    \item
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{https://www.csie.ntu.edu.tw/~cjlin/libsvm/}

    \item Ming-Wei Chang, Hsuan-Tien Lin, Ming-Hen Tsai, Chia-Hua Ho and Hsiang-Fu Yu\cr
     \emph{Weights for data instances}\cr
     \url{https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/#weights_for_data_instances}

     \item
      Exact formulations of models, algorithms, etc. can be found in the
      document:\cr
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{https://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.ps.gz}

      \item
      More implementation details and speed benchmarks can be found on:
      Rong-En Fan and Pai-Hsune Chen and Chih-Jen Lin:\cr
      \emph{Working Set Selection Using the Second Order Information for Training SVM}\cr
      \url{https://www.csie.ntu.edu.tw/~cjlin/papers/quadworkset.pdf}

  }
}
\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{
  \code{\link{predict.wsvm}},
  \code{\link{plot.wsvm}},
  \code{\link{tune_wsvm}},
  \code{\link[SparseM]{matrix.csr}} (in package \pkg{SparseM})
}
\examples{
d
}










