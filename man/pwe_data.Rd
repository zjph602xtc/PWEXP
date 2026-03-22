\name{pwe_data}
\alias{pwe_data}
\title{Piece-wise Exponentail Example Data}
\usage{pwe_data
       data(pwe_data, package="PwePred")
}
\description{This is a simulated clinical data with piece-wise exponential distributed endpoint.}
\format{
\tabular{ll}{
ID:\tab id\cr
randT:\tab the randomization time for each subject\cr
eventT:\tab the event time\cr
dropT:\tab the drop-out time\cr
censor_reason:\tab If a subject is censored, \code{censor_reason} shows the type of censoring (i.e., \code{drop_out}, \code{death})\cr
event:\tab indicates whether the primary event occurred at the end of follow-up with $0$ for censoring and $1$ for occurrence of event\cr
followT:\tab the follow-up time, which is the minimum value of \code{eventT}, \code{dropT}. In real-world datasets, this is the observation time\cr
followT_abs:\tab the sum of \code{randT} and \code{followT}
}}
\note{Here we utilize the \code{simdata} function to create a simple example dataset with the following characteristics:
\itemize{
  \item randomization rate is defined as 20 subjects per month (\code{rand_rate = 20}), and total sample size is 1000 (\code{total_sample = 1000});
  \item the primary endpoint (event)  follows a PWE distribution with monthly hazard rates of 0.1 for t < 5 months, 0.01 for 5 <= t < 14 months, 0.2 for t >= 14 months (defined by \code{myevent_dist} function);
  \item the drop-out follows an exponential distribution with a drop-out probability of 3\%/month (\code{drop_rate = 0.03}) (equivalently, a monthly drop-out hazard rate is -log(1-0.03)=0.0304);
  \item the \code{add_column} argument requests additional variables to be included in the dataset.
}
\code{set.seed(1818)}

\code{myevent_dist <- function(n)rpwexpm(n, rate = c(0.1, 0.01, 0.2), breakpoint = c(5,14))}

\code{pwe_data <- simdata(rand_rate = 20, total_sample = 1000, drop_rate = 0.03, advanced_dist = list(event_dist = myevent_dist), add_column = c('censor_reason', 'event', 'followT', 'followT_abs'))}
}


\keyword{PwePred}
