plot.pwexp.fit <- function(x, ...){
  stop('Please use \'plot_survival\' function to visualize the model.')
}

plot.boot.pwexp.fit <- function(x, ...){
  stop('Please use \'plot_survival\' function to visualize the model.')
}

plot_survival <- function (time, ...){
  UseMethod ("plot_survival")
}

plot_survival.default <- function(time, event, add=FALSE, conf.int=FALSE, mark.time=TRUE,
                                  lwd=2, xlab='Follow-up time', ylab='Survival function', ...){
  arg <- list(...)
  option <- c('conf.int', 'mark.time', 'lwd', 'xlab', 'ylab')
  default <- list(conf.int, mark.time, lwd, xlab, ylab)
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]

  if (is.data.frame(time)){
    time <- time[[1]]
  }
  if (is.data.frame(event)){
    event <- event[[1]]
  }

  s <- survfit(Surv(time, event)~1)
  do.call(ifelse(add, 'lines', 'plot'), c(list(x=s), arg))
}


plot_survival.pwexp.fit <- function(time, add=TRUE, show_breakpoint=TRUE,
                                    breakpoint_par=NULL, ...){
  object <- time
  arg <- list(...)
  option <- c('lwd', 'col', 'xlab', 'ylab','type')
  default <- list(2, 'red', 'Follow-up time', 'Survival function','l')
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]

  option_brk <- c('lwd', 'col', 'lty')
  default_brk <- list(2, 'grey', 2)
  ind_brk <- option_brk %in% names(breakpoint_par)
  breakpoint_par[option_brk[!ind_brk]] <- default_brk[!ind_brk]

  if (add){
    xrange <- seq(0, par('usr')[2], length=200)
    y <- PWEXP::ppwexp(xrange, rate=attr(object,'lam'),
                breakpoint=attr(object,'brk'), lower.tail = F)
    do.call(lines, c(list(x=xrange, y=y), arg))
  }else{
    if (!is.null(arg$xlim)){
      xrange <- seq(0, arg$xlim[2], length=200)
    }else if (!is.null(attr(object,'brk'))){
      xrange <- seq(0, max(attr(object,'brk'))*1.4, length=200)
    }else {
      xrange <- seq(0, 2/min(attr(object,'lam')), length=200)
    }
    y <- PWEXP::ppwexp(xrange, rate=attr(object,'lam'),
                breakpoint=attr(object,'brk'), lower.tail = F)
    do.call(plot, c(list(x=xrange, y=y), arg))
  }
  if (show_breakpoint && !is.null(attr(object,'brk'))){
    do.call(abline, c(list(v=attr(object,'brk')), breakpoint_par))
  }
}

plot_survival.boot.pwexp.fit <- function(time, add=TRUE, alpha=0.1, show_breakpoint=TRUE,
                                         breakpoint_par=NULL, show_CI=TRUE, CI_par=NULL, ...){
  # arg <- list(...)
  # option <- c('lwd', 'xlab', 'ylab')
  # default <- list(2, 'Follow-up time', 'Survival function')
  # ind <- option %in% names(arg)
  # arg[option[!ind]] <- default[!ind]
  #
  # option_brk <- c('lwd', 'col', 'lty')
  # default_brk <- list(2, 'grey', 2)
  # ind_brk <- option_brk %in% names(breakpoint_par)
  # breakpoint_par[option_brk[!ind_brk]] <- default_brk[!ind_brk]

  obj <- object <- time
  obj <- obj[1,,drop=F]
  attr(obj, 'lam') <- as.numeric(attr(obj, 'lam')[1,])
  if (!is.null(attr(obj, 'brk'))){
    attr(obj, 'brk') <- as.numeric(attr(obj, 'brk')[1,])
  }
  class(obj) <- c('pwexp.fit','data.frame')
  plot_survival.pwexp.fit(obj, add=add, show_breakpoint = show_breakpoint, breakpoint_par = breakpoint_par, ...)

  option_ci <- c('lwd', 'col', 'lty')
  default_ci <- list(2, '#ff9896', 2)
  ind_ci <- option_ci %in% names(CI_par)
  CI_par[option_ci[!ind_ci]] <- default_ci[!ind_ci]

  if (show_CI){
    xrange <- seq(0, par('usr')[2], length=200)
    if (!is.null(attr(object,'brk'))){
      line_data <- mapply(function(rate, breakpoint)PWEXP::ppwexp(xrange, rate = rate, breakpoint = breakpoint, lower.tail = F),
                          rate=as.list(as.data.frame(t(attr(object,'lam')))),
                          breakpoint=as.list(as.data.frame(t(attr(object,'brk')))))
    }else{
      line_data <- mapply(function(rate)PWEXP::ppwexp(xrange, rate = rate, lower.tail = F),
                          rate=as.list(as.data.frame(t(attr(object,'lam')))))
    }
    ci_data <- apply(line_data, 1, function(x)quantile(x, c(alpha/2, (1-alpha/2)), na.rm=T))
    do.call(lines, c(list(x=xrange, y=ci_data[1,]), CI_par))
    do.call(lines, c(list(x=xrange, y=ci_data[2,]), CI_par))
  }
}

plot_event <- function (time, ...){
  UseMethod ("plot_event")
}

plot_event.default <- function(time, event, abs_time=TRUE, additional_event=0, add=FALSE,
                               plot=TRUE, xyswitch=FALSE, ...){
  # here dat$followT_abs may be followT or followT_abs, depending on the user input
  arg <- list(...)
  if (xyswitch){
    option <- c('lwd', 'xlab', 'ylab', 'type')
    default <- list(2,'Cumulative events', ifelse(abs_time, 'Time from first randomization', 'Time from randomization of each subject'), 'l')
  }else{
    option <- c('lwd', 'xlab', 'ylab', 'type')
    default <- list(2, ifelse(abs_time, 'Time from first randomization', 'Time from randomization of each subject'), 'Cumulative events', 'l')
  }
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]

  if (is.data.frame(time)){
    time <- time[[1]]
  }
  if (is.data.frame(event)){
    event <- event[[1]]
  }

  dat <- data.frame(followT_abs = time, event = event)
  dat <- dat[order(dat$followT_abs),]
  dat$cum <- cumsum(dat$event)+additional_event
  if (plot){
    if (xyswitch){
      do.call(ifelse(add, 'lines', 'plot'), c(list(x=dat$cum, y=dat$followT_abs), arg))
    }else{
      do.call(ifelse(add, 'lines', 'plot'), c(list(x=dat$followT_abs, y=dat$cum), arg))
    }
  }else{
    line_data <- data.frame(time=dat$followT_abs, n_event=dat$cum)
    return(line_data)
  }
}



plot_event.predict.pwexp.fit <- function(time, abs_time=TRUE, add=TRUE, plot=TRUE, xyswitch=FALSE,
                                         eval_at=NULL, ...){
  predict_model <- time
  arg <- list(...)
  if (xyswitch){
    option <- c('lwd', 'col', 'xlab', 'ylab', 'type')
    default <- list(2, 'red', 'Cumulative events', 'Time from first randomization', 'l')
  }else{
    option <- c('lwd', 'col', 'xlab', 'ylab', 'type')
    default <- list(2, 'red', 'Time from first randomization', 'Cumulative events', 'l')
  }
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]

  if (!abs_time){
    stop('abs_time must be TRUE when plotting the prediction curve. ')
  }

  #  to predict type='event' or 'time'
  if (plot){
    if (add){
      xrange <- seq(0, par('usr')[2], length=200)
    }else{
      if (!is.null(arg$xlim)){
        xrange <- seq(0, arg$xlim[2], length=200)
      }else {
        stop('Must provide xlim when NOT adding the prediction curve to an existing figure. ')
      }
    }
    x_pre_range <- seq(0, max(xrange[200], 240), length=5000)
  }else{
    x_pre_range <- seq(0, 240, length=10000)
  }
  pre <- sapply(predict_model$event_fun, function(f)f(x_pre_range))
  if (!is.matrix(pre) | is.data.frame(pre)){
    pre <- matrix(pre, nrow=1)
  }
  # na_include <- rowMeans(is.na(pre)) < 0.05
  pre <- apply(pre, 1, function(x)quantile(x, 0.5,na.rm=T))
  # pre[,!na_include] <- NA

  if (!xyswitch){
    median_line <- suppressWarnings(approxfun(x_pre_range, pre, rule=1:2, ties='min'))
  }else{
    median_line <- function(x){
      res <- rep(NA, length=length(x))
      tmp_ind <- x < max(pre, na.rm = T)
      res[tmp_ind] <- suppressWarnings(approxfun(pre, x_pre_range, rule=1, ties='min')(x[tmp_ind]))
      return(res)
    }
  }

  if (plot){
    do.call(ifelse(add, 'lines', 'plot'), c(list(x=xrange, y=median_line(xrange)), arg))
  }

  if (!is.null(eval_at)){
    pre <- cbind(eval_at, median_line(eval_at))
    if (xyswitch){
      colnames(pre) <- c('n_event', 'time')
    }else{
      colnames(pre) <- c('time','n_event')
    }
    return(pre)
  }
}

plot_event.predict.boot.pwexp.fit <- function(time, abs_time=TRUE,  alpha=0.1, type='confidence', add=TRUE, plot=TRUE, xyswitch=FALSE, eval_at=NULL, show_CI=TRUE, CI_par=NULL, ...){
  predict_model <- time
  arg <- list(...)
  if (xyswitch){
    option <- c('lwd', 'col', 'xlab', 'ylab', 'type')
    default <- list(2, 'red', 'Cumulative events', 'Time from first randomization', 'l')
  }else{
    option <- c('lwd', 'col', 'xlab', 'ylab', 'type')
    default <- list(2, 'red', 'Time from first randomization', 'Cumulative events', 'l')
  }
  ind <- option %in% names(arg)
  arg[option[!ind]] <- default[!ind]

  option_ci <- c('lwd', 'col', 'lty')
  default_ci <- list(2, '#ff9896', 2)
  ind_ci <- option_ci %in% names(CI_par)
  CI_par[option_ci[!ind_ci]] <- default_ci[!ind_ci]

  if (!abs_time){
    stop('abs_time must be TRUE when plotting the prediction curve. ')
  }
  if (!(type %in% c('predictive','confidence'))){
    stop('wrong type argument.')
  }

  #  to predict type='event' or 'time'
  if (plot){
    if (add){
      xrange <- seq(0, par('usr')[2], length=200)
    }else{
      if (!is.null(arg$xlim)){
        xrange <- seq(0, arg$xlim[2], length=200)
      }else {
        stop('Must provide xlim when NOT adding the prediction curve to an existing figure. ')
      }
    }
    x_pre_range <- seq(0, max(xrange[200], 240), length=1000)
  }else{
    x_pre_range <- seq(0, 240, length=2000)
  }
  pre <- sapply(predict_model$event_fun, function(f)f(x_pre_range))
  if (!is.matrix(pre) | is.data.frame(pre)){
    pre <- matrix(pre, ncol=1)
  }
  # na_include <- rowMeans(is.na(pre)) < 0.05
  if (type=='confidence'){
    grp_con <- sapply(strsplit(colnames(pre), split = '_'), '[', 3)
    pre <- aggregate(as.data.frame(t(pre)), by = list(grp_con), function(x)mean(x, na.rm = T), drop=F)[,-1,drop=F]
  }else{
    pre <- t(pre)
  }
  pre <- apply(pre, 2, function(x)quantile(x, c(alpha/2, 0.5, (1-alpha/2)),na.rm=T))
  # pre[,!na_include] <- NA

  if (!xyswitch){
    median_line <- suppressWarnings(approxfun(x_pre_range, pre[2,], rule=1:2, ties='min'))
    low_line <- suppressWarnings(approxfun(x_pre_range, pre[1,], rule=1:2, ties='min'))
    up_line <- suppressWarnings(approxfun(x_pre_range, pre[3,], rule=1:2, ties='min'))
  }else{
    median_line <- function(x){
      res <- rep(NA, length=length(x))
      tmp_ind <- x < max(pre[2,], na.rm = T)
      res[tmp_ind] <- suppressWarnings(approxfun(pre[2,], x_pre_range, rule=1, ties='min')(x[tmp_ind]))
      return(res)
    }
    low_line <- function(x){
      res <- rep(NA, length=length(x))
      tmp_ind <- x < max(pre[3,], na.rm = T)
      res[tmp_ind] <- suppressWarnings(approxfun(pre[3,], x_pre_range, rule=1, ties='min')(x[tmp_ind]))
      return(res)
    }
    up_line <- function(x){
      res <- rep(NA, length=length(x))
      tmp_ind <- x < max(pre[1,], na.rm = T)
      res[tmp_ind] <- suppressWarnings(approxfun(pre[1,], x_pre_range, rule=1, ties='min')(x[tmp_ind]))
      return(res)
    }
  }

  if (plot){
    do.call(ifelse(add, 'lines', 'plot'), c(list(x=xrange, y=median_line(xrange)), arg))
    if (show_CI){
      do.call(lines, c(list(x=xrange, y=low_line(xrange)), CI_par))
      do.call(lines, c(list(x=xrange, y=up_line(xrange)), CI_par))
    }
  }

  if (!is.null(eval_at)){
    pre <- cbind(eval_at, median_line(eval_at), low_line(eval_at), up_line(eval_at))
    if (xyswitch){
      colnames(pre) <- c('n_event', 'time', paste0(alpha/2*100, '% time'), paste0((1-alpha/2)*100, '% time'))
    }else{
      colnames(pre) <- c('time','n_event', paste0(alpha/2*100, '% n_event'), paste0((1-alpha/2)*100, '% n_event'))
    }
    return(pre)
  }

}
