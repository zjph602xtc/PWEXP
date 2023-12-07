\name{NEWS}
\title{News for Package \pkg{PWEXP}}

\section{Changes in version 0.4.2}{
  \itemize{
    \item(pwexp.fit) add empirical number of breakpoints when 'nbreak=NULL'
    \item(vignettes) add 'Usage of Package'
    \item(pwexp.fit, boot.pwexp.fit, cv.pwexp.fit) add 'min_pt_tail' parameter
    \item(sim_followup) change parameter name 'min_follow' to 'extra_follow'
    \item(sim_followup) change parameter name 'count_in_min_follow' to 'count_in_extra_follow'
  }
}

\section{Changes in version 0.4.1}{
  \itemize{
    \item (plot_event) rewrite 'xyswitch=TRUE' to have the same result as 'xyswitch=FALSE'
    \item add 'Depends: survival, fastmatch' to DESCRIPTION
    \item (boot.pwexp.fit, cv.pwexp.fit) correctly attach functions in 'doSNOW', 'foreach' packages
  }
}

\section{Changes in version 0.4.0}{
  \itemize{
    \item (boot.pwexp.fit) add 'parallel' option
    \item (cv.pwexp.fit) add 'parallel' option
  }
}

\section{Changes in version 0.3.1}{
  \itemize{
    \item (plot_event) fix an error that 'xlab', 'ylab' do not work when 'xyswitch=TRUE'
  }
}

\section{Changes in version 0.3.0}{
  \itemize{
    \item (sim_followup) fix an error when type='sample' or 'event', but 'at' is not an integer
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

\section{Changes in version 0.1.0}{
  \itemize{
    \item The first release.
  }
}
