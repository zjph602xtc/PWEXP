\name{NEWS}
\title{News for Package \pkg{PWEXP}}

\section{Changes in version 0.1.0}{
  \itemize{
    \item The first release.
  }
}

\section{Changes in version 0.2.0}{
  \itemize{
    \item (pwexp.fit) fix when number of breaks is 2 more larger than number of event (error in combn function)
    \item (pwexp.fit) add warning 'Incorrect result returned'
    \item (pwexp.fit) fix case when no events and nbreak=0, now return -Inf results instead of 0
    \item (pwexp.fit) add 'tol' parameter to control the minimum gap between two breakpoints
    \item (predict.boot.pwexp.fit) skip unsuccessful prediction
    \item (predict.boot.pwexp.fit) add first row c(analysis_time, sum(para$event)) to line_data to make prediction from analysis_time
    \item (sim_followup) fix number of event=0 error
    \item (sim_followup) fix number of enrolled=0 error
  }
}

\section{Changes in version 0.3.0}{
  \itemize{
    \item (sim_followup) fix an error when type='sample' or 'event', but 'at' is not an integer
  }
}
