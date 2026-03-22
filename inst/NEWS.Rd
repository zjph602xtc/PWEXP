\name{NEWS}
\title{News for Package \pkg{PwePred}}

\section{Changes in version 1.1.1}{
  \itemize{
    \item (AIC.pwexpm, BIC.pwexpm, logLik.pwexpm, confint.boot.pwexpm) export these functions in the namespace
    \item (pew_data, pwe_data_robust_tail, ruppert_ipd) add example datasets
  }
}

\section{Changes in version 1.1.0}{
  \itemize{
    \item (summary.pwexpm, print.summary.pwexpm) print the actual optimizer used for estimation
    \item (summary.boot.pwexpm, print.summary.boot.pwexpm) print the actual optimizer used for estimation; print the number of successful iterations
    \item (summary.predict.pwexpm, print.summary.predict.pwexpm, summary.predict.boot.pwexpm, print.summary.predict.boot.pwexpm) print the number of successful simulations
    \item (coef.pwexpm) extract coefficents in 'pwexpm'
    \item (confint.boot.pwexpm) confidence intervals for model parameters
    \item (sim_followup) update the help document; add 'sim_followup' class
    \item (print.sim_followup) add the print function for 'sim_followup' class
  }
}

\section{Changes in version 1.0.0}{
  \itemize{
    \item name the package as 'PwePred'
    \item (pwexpm, boot.pwexpm, cv.pwexpm) add 'pwexpm','boot.pwexpm','cv.pwexpm' functions that use Surv/data as inputs
    \item change the object name 'pwexpm.fit' to 'pwexpm', and removed all '.fit' part from all other objects.
    \item (dpwexpm, ppwexpm, qpwexpm, rpwexpm) change the function names from 'pwexp', 'ppwexp', 'qpwexp', 'rpwexp' to avoid replicate with the 'rpact' package

    \item (pwexpm_fit, pwexpm) return a list instead of a matrix
    \item (print.pwexpm, AIC.pwexpm, BIC.pwexpm, logLiK.pwexpm) add S3 methods for 'pwexpm' object
    \item (plot.pwexpm) add default plot for 'pwexpm' object
    \item (plot_survival.pwexpm) return stop if no previous plot is called

    \item (boot.pwexpm_fit, boot.pwexpm) return a list instead of a matrix
    \item (print.boot.pwexpm) add S3 methods for 'boot.pwexpm' object
    \item (plot.boot.pwexpm) add default plot for 'boot.pwexpm' object
    \item (plot_survival.boot.pwexpm) return stop if no previous plot is called

    \item (boot.pwexpm.pwexpm, cv.pwexpm.pwexpm) update to fit the new 'pwexpm' object

    \item (plot.cv.pwexpm) add default plot for 'cv.pwexpm' object
    \item (print.cv.pwexpm) add S3 methods for 'print.cv.pwexpm' object

    \item (predict.pwexpm, predict.boot.pwexpm) update to fit the new 'pwexpm', 'boot.pwexpm' functions
    \item (predict.pwexpm, predict.boot.pwexpm) return a list with more model details
    \item (print.predict.pwexpm, print.predict.boot.pwexpm) add S3 methods for 'predict.pwexpm', 'predict.boot.pwexpm' object
    \item (predict.boot.pwexpm) change parameter name 'censor_model_boot' to 'censor_model'
    \item (plot.predict.boot.pwexpm, plot.predict.pwexpm) add default plot
    \item (simdata) set add_column=c('followT','event') as the default value

    \item (pwexpm, pwexpm_fit, cv.pwexpm, cv.pwexpm_fit, boot.pwexpm, boot.pwexpm_fit, sim_followup, predict.pwexpm, predict.boot.pwexpm) allow 'seed=NULL'
  }
}

\section{Changes in version 0.5.0}{
  \itemize{
    \item(predict.boot.pwexp.fit) remove 'type' parameter
    \item(plot_event.predict.boot.pwexp.fit) add 'type' parameter; change the dafault to 'confidence'
  }
}

\section{Changes in version 0.4.5}{
  \itemize{
    \item(sim_followup) add 'follow_up_endpoint' parameter, which is used to define the end of follow-up
  }
}

\section{Changes in version 0.4.4}{
  \itemize{
    \item(simdata) allow non-integer values in 'n_rand', 'rand_rate'
    \item(sim_followup) allow non-integer values in 'n_rand', 'rand_rate'
  }
}

\section{Changes in version 0.4.3}{
  \itemize{
    \item(simdata) change enroll behavior. For example, when rand_rate=20, total_sample=25, in the previous version, randT is between 0 to 2, but now randT is between 0 to 1.25.
  }
}

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
