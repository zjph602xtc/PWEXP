plot_km <- function(time, event, ...){
  arg <- list(...)
  option <- c('conf.int', 'mark.time', 'lwd', 'xlab', 'ylab')
  default <- list(FALSE, TRUE, 2, 'Follow-up time', 'Survival function')
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]
  
  s <- survfit(Surv(time, event)~1)
  do.call(plot, c(list(x=s), arg))
}


plot_event <- function(time, abs_time=TRUE, event, additional_event=0, add=FALSE, 
                       plot=TRUE, xyswitch=FALSE, ...){
  # here dat$followT_abs may be followT or followT_abs, depending on the user input
  arg <- list(...)
  option <- c('lwd', 'xlab', 'ylab', 'type')
  default <- list(2, ifelse(abs_time, 'Time from the start of the trial', 'Individual time 
                            from enrollment'), 'Cumulative events', 'l')
  if (xyswitch){
    arg$ylab <- arg$xlab
    arg$xlab <- 'Cumulative events'
  }
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]
  
  dat <- data.frame(followT_abs = time, event = event)
  dat <- dat[order(dat$followT_abs),]
  dat$cum <- cumsum(dat$event)+additional_event
  if (plot){
    if (xyswitch){
      do.call(ifelse(add, 'lines', 'plot'), c(list(x=dat$cum, y=dat$followT_abs), arg))
    }else{
      do.call(ifelse(add, 'lines', 'plot'), c(list(x=dat$followT_abs, y=dat$cum), arg))
    }
  }
  
  line_data <- data.frame(time=dat$followT_abs, n_event=dat$cum)
  return(line_data)
}
