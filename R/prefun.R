library(survival)
library(fastmatch)

rpwexp <- function(n, rate, breakpoint=NULL){
  x <- runif(n)
  return(qpwexp(x, rate, breakpoint, lower.tail=TRUE, log.p=FALSE))
}

rpwexp_conditional <- function(n, qT, rate, breakpoint=NULL){
  # length of qT should be 1 or == length(n)
  x <- runif(n)
  return(qpwexp_conditional(x, qT, rate, breakpoint, lower.tail=TRUE, log.p=FALSE))
}

ppwexp <- function(q, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # S(t)=exp(-((lam1-lam2)*d1+(lam2-lam3)*d2+... lami*t))
  # breakpoint must be sorted!
  # if safety_check is FALSE, then skip all checks and one_piece must be correctly specified. 
  if (safety_check){
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
  }
  if (one_piece){
    return(pexp(q, rate, lower.tail, log.p))
  }else {
    shift_t <- c(0, cumsum(-diff(rate)*breakpoint))
    interval <- findInterval(q, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
    logs <- -(rate[interval]*q+shift_t[interval])
    if (!lower.tail & log.p){
      return(logs)
    }else if (!lower.tail & !log.p){
      return(exp(logs))
    }else if (lower.tail & !log.p){
      return(1-exp(logs))
    }else {
      return(log(1-exp(logs)))
    }
  }
}

ppwexp_conditional <- function(q, qT, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # F(t|X>qT) = 1 - S(t)/S(qT)
  # S(t|X>qT) = S(t)/S(qT)
  # breakpoint must be sorted!
  # if safety_check is FALSE, then skip all checks and one_piece must be correctly specified. 
  if (safety_check){
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
    if (length(qT)!=1 & length(q)!=length(qT)){
      stop('length of qT should be 1 or == length(q)')
    }
    if (any(q < qT)){
      stop('q must be equal or larger than qT')
    }
  }
  logs <- ppwexp(q, rate, breakpoint, lower.tail = F, log.p = T, one_piece, safety_check = FALSE)-
    ppwexp(qT, rate, breakpoint, lower.tail = F, log.p = T, one_piece, safety_check = FALSE)
  if (!lower.tail & log.p){
    return(logs)
  }else if (!lower.tail & !log.p){
    return(exp(logs))
  }else if (lower.tail & !log.p){
    return(1-exp(logs))
  }else {
    return(log(1-exp(logs)))
  }
}

dpwexp <- function(x, rate=1, breakpoint=NULL, log=FALSE, one_piece, safety_check=TRUE){
  # f(t)=lambda_i * exp(-((lam1-lam2)*d1+(lam2-lam3)*d2+... lami*t))
  # breakpoint must be sorted!
  if (safety_check){
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
  }
  if (one_piece){
    return(dexp(x, rate, log))
  }else {
    shift_t <- c(0, cumsum(-diff(rate)*breakpoint))
    log_rate <- log(rate)
    interval <- findInterval(x, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
    logd <- -rate[interval]*x-(shift_t-log_rate)[interval]
    # logd <- -(rate[interval]*x+shift_t[interval])+log(rate[interval])
    if (log){
      return(logd)
    }else {
      return(exp(logd))
    }
  }
}

qpwexp <- function(p, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # F(y)^(-1)=(-log(1-y)-(lam1-lam2)*d1-(lam2-lam3)*d2-...)/lami
  # breakpoint must be sorted!
  if (safety_check){
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
  }
  if (one_piece){
    return(qexp(p, rate, lower.tail, log.p))
  }else {
    if (log.p){
      p <- exp(p)
    }
    shift_t <- c(0, cumsum(-diff(rate)*breakpoint))
    Fcut <- ppwexp(breakpoint, rate, breakpoint, lower.tail, log.p=FALSE, one_piece, safety_check = FALSE)
    if (lower.tail){
      interval <- findInterval(p, vec=c(-Inf, Fcut, Inf), rightmost.closed = FALSE)
      t <- (-log(1-p)-shift_t[interval])/rate[interval]
    }else{
      interval <- 2+length(Fcut)-findInterval(p, vec=rev(c(Inf, Fcut, -Inf)), left.open = TRUE, rightmost.closed = TRUE)
      t <- (-log(p)-shift_t[interval])/rate[interval]
    }
    return(t)
  }
}

qpwexp_conditional <- function(p, qT, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # t=F(y)^(-1) given X>=qT
  # assume qT in Mth interval, t (or say p) in kth interval
  # t=(-log(1-y)-(lam_{M}-lam_{M+1})*d_M-...-(lam_{k-1}-lam_{k})*d_{k-1}+lam_M*T)/lam_{K}
  # breakpoint must be sorted!
  
  # length of qT should be 1 or == length(p)
  if (safety_check){
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
    if (length(qT)!=1 & length(p)!=length(qT)){
      stop('length of qT should be 1 or == length(p)')
    }
  }
  if (log.p){
    p <- exp(p)
  }
  if (one_piece){
    if (lower.tail){
      t <- (-log(1-p)+rate*qT)/rate
    }else {
      t <- (-log(p)+rate*qT)/rate
    }
    return(t)
  }
  if (length(qT)!=1){
    res <- mapply(function(p,qT)qpwexp_conditional(p, qT, rate, breakpoint, lower.tail, log.p, one_piece, safety_check = FALSE), p=p,qT=qT)
    return(res)
  }else {
    shift_t <- -diff(rate)*breakpoint
    Fcut <- ppwexp_conditional(breakpoint, qT, rate, breakpoint, lower.tail, log.p=FALSE, one_piece, safety_check = FALSE)
    if (lower.tail){
      interval_k <- findInterval(p, vec=c(-Inf, Fcut, Inf), rightmost.closed = FALSE)
      interval_m <- findInterval(qT, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
      # shift_int <- mapply(function(i,j)shift_t[i:j], i=interval_m, j=interval_k-1, SIMPLIFY = FALSE)
      shift_int <- lapply(interval_k-1, function(k)shift_t[interval_m:k])
      shift_int[interval_m > (interval_k-1)] <- 0
      shift_int <- sapply(shift_int, sum)
      t <- (-log(1-p)-shift_int+rate[interval_m]*qT)/rate[interval_k]
    }else{
      interval_k <- 2+length(Fcut)-findInterval(p, vec=rev(c(Inf, Fcut, -Inf)), left.open = TRUE, rightmost.closed = TRUE)
      interval_m <- findInterval(qT, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
      # shift_int <- mapply(function(i,j)shift_t[i:j], i=interval_m, j=interval_k-1, SIMPLIFY = FALSE)
      shift_int <- lapply(interval_k-1, function(k)shift_t[interval_m:k])
      shift_int[interval_m > (interval_k-1)] <- 0
      shift_int <- sapply(shift_int, sum)
      t <- (-log(p)-shift_int+rate[interval_m]*qT)/rate[interval_k]
    }
    return(t)  
  }
}

get_grid <- function(time, breakpoint, nbreak, max_set=5000){
  time <- unique(time)
  time <- time[-1]
  
  if (nbreak==0){
    nbreak <- length(breakpoint)
  }
  if (!is.null(breakpoint)){
    nbreak <- nbreak - length(breakpoint)
    if (nbreak==0){
      return(breakpoint <- matrix(breakpoint, nrow=1))
    }
    if (nbreak<0){
      stop('nbreak is the total number of breakpoints, which must be equal or larger than the length of \'breakpoint\'')
    }
    time <- setdiff(time, breakpoint)
    except <- sapply(breakpoint, function(x) findInterval(x, time))
    time <- time[-c(except, except+1)]
  }
  nr <- length(time)
  nl <- 1
  
  if (choose(nr, nbreak) > max_set){
    while((nr-nl)>1.1){
      if (choose(floor((nl+nr)/2), nbreak) > max_set){
        nr <- floor((nl+nr)/2)
      }else{
        nl <- floor((nl+nr)/2)
      }
    }
  }
  time <- sort(sample(time,nr,replace = F))
  setting <- t(combn(time, nbreak))
  
  if (NROW(setting) > max_set){
    setting <- setting[sort(sample.int(NROW(setting), max_set, replace = F)),,drop=F]
  }
  
  if (!is.null(breakpoint)){
    setting <- t(apply(setting, 1, function(x)sort(c(breakpoint, x))))
  }
  return(setting)
}

fit_pwexp <- function(time, event, breakpoint=NULL, nbreak=0, max_set=10000, seed=1818, trace=FALSE){
  # n_{fj}/lam_j - sum_{set fj sj}(t_i-d_{j-1}) = n_{set f_j+1 to f_end, s_j+1 to s_end}*(d_j-d_{j-1})
  # f is event data, s is censored data
  
  # auto search if breakpoint=NULL
  ind <- order(time)
  event <- event[ind]
  time <- time[ind]
  time_event <- time[event==1]
  time_noevent <- time[event==0]
  N <- length(time)
  set.seed(seed)
  n_fix_brk <- length(breakpoint)
  if (is.null(breakpoint)  & nbreak==0){
    lam <- sum(event)/sum(time)
    loglikelihood <- sum(dpwexp(time_event, rate=lam, breakpoint = NULL, log = T, one_piece = T, safety_check = F))+
      sum(ppwexp(time_noevent, rate=lam, lower.tail = F, breakpoint = NULL, log.p = T, one_piece = T, safety_check = F))
    aic <- 2-2*loglikelihood
    bic <- log(length(time))-2*loglikelihood
    res <- data.frame(lam1=lam, likelihood=loglikelihood, AIC=aic, BIC=bic)
    attr(res,'lam') <- lam
    return(res)
  }
  if (!is.null(breakpoint) & max(breakpoint) > max(time_event)){
      warning('One of the \'breakpoint\' is too large. No event after this breakpoint. We will NOT use it. ')
      breakpoint <- breakpoint[breakpoint <= max(time_event)]
      n_fix_brk <- length(breakpoint)
      if (n_fix_brk==0){
        breakpoint <- NULL
      }
  }

  breakpoint <- get_grid(time, breakpoint, nbreak, max_set)
  
  res <- matrix(-Inf, ncol=2*NCOL(breakpoint)+2, nrow = NROW(breakpoint))
  for (i in 1:NROW(breakpoint)){
    brk0 <- breakpoint[i,]
    brk <- c(0, brk0, Inf)
    tmpi <- findInterval(time, vec=brk)
    numerator <- ctapply(event, tmpi, sum)
    if (any(numerator==0)){
      next
    }
    # rhs <- c(diff(c(0, brk0)) * sapply(brk0, function(x)sum(time >= x)),0)
    rhs <- c(diff(c(0, brk0)),0) * (N-cumsum(ctapply(time, tmpi, length)))
    # lhs_term2 <- sapply(1:(length(brk)-1), function(x)sum(time[time < brk[x+1] & time >= brk[x]] - brk[x]))
    lhs_term2 <- ctapply(time - brk[tmpi], tmpi, sum)
    # lam <- sapply(1:(length(brk)-1), function(x)sum(time < brk[x+1] & time >= brk[x] & event))/(rhs+lhs_term2)
    lam <- numerator/(rhs+lhs_term2)
    loglikelihood <- sum(dpwexp(time_event, rate=lam, breakpoint = brk0, log = T, one_piece = F, safety_check = F))+
      sum(ppwexp(time_noevent, rate=lam, lower.tail = F, breakpoint = brk0, log.p = T, one_piece = F, safety_check = F))
    res[i,] <- c(brk0, lam, loglikelihood)
  }
  res <- data.frame(res)
  if (!trace){
    res <- res[which.max(res[,NCOL(res)])[1],,drop=F]
    n_k <- 2*max(n_fix_brk, nbreak)+1
    aic <- 2*(n_k-n_fix_brk)-2*res[,NCOL(res)]
    bic <- (n_k-n_fix_brk)*log(N)-2*res[,NCOL(res)]
    res <- cbind(res, aic, bic)
    attr(res,'lam') <- as.numeric(res[,(NCOL(breakpoint)+1):(2*NCOL(breakpoint)+1)])
    attr(res,'brk') <- as.numeric(res[,1:NCOL(breakpoint)])
  }
  colnames(res) <- c(paste0('brk', 1:NCOL(breakpoint)), paste0('lam', 1:(NCOL(breakpoint)+1)), 'likelihood','AIC','BIC')
  return(res)
  
}

plot_km <- function(dat, ...){
  plot(survfit(Surv(followT, event)~1, data=dat), conf.int = F, mark.time = T, lwd=2, xlab='abs time from start', ...)
}

plot_event <- function(dat, add=F, abs_time=T, additional_event=0, plot=TRUE, xyswitch=F, ...){
  if (abs_time){
    dat <- dat[order(dat$followT_abs),]
    dat$cum <- cumsum(dat$event)+additional_event
    if (plot){
      if (add){
        if (xyswitch){
          lines(dat$cum, dat$followT_abs, lwd=2, ...)
        }else{
          lines(dat$followT_abs, dat$cum, lwd=2, ...)
        }
      }else{
        if (xyswitch){
          plot(dat$cum, dat$followT_abs, type='l',lwd=2, ...)
        }else {
          plot(dat$followT_abs, dat$cum, type='l',lwd=2, ...)
        }
      }
    }
    line_data <- data.frame(followT_abs=dat$followT_abs, n_event=dat$cum)
  }else{
    dat <- dat[order(dat$followT),]
    dat$cum <- cumsum(dat$event)+additional_event
    if (plot){
      if (add){
        if (xyswitch){
          lines(dat$cum, dat$followT, lwd=2, ...)
        }else {
          lines(dat$followT, dat$cum, lwd=2, ...)
        }
      }else{
        if (xyswitch){
          plot(dat$cum, dat$followT, type='l',lwd=2, ...)
        }else {
          plot(dat$followT, dat$cum, type='l',lwd=2, ...)
        }
      }
    }
    line_data <- data.frame(followT=dat$followT, n_event=dat$cum)
  }
  return(line_data)
}

cut_dat <- function(dat, cut){
  train <- dat[dat$enrollT <= cut,]
  train$censor_reason[train$followT_abs > cut] <- 'cut'
  train$event[train$followT_abs > cut] <- 0
  train$followT_abs[train$followT_abs > cut] <- cut
  train$followT <- train$followT_abs-train$enrollT
  return(train)
}

tran_censor <- function(dat){
  dat$event <- ifelse(dat$censor_reason=='censored', 1, 0)
  dat$event[is.na(dat$event)] <- 0
  return(dat)
}

fit_pwexp_bootstrap <- function(time, event, nsim=100, breakpoint=NULL, nbreak=0, max_set=1000, seed=1818){
  dat <- data.frame(time=time, event=event)
  n <- NROW(dat)
  res_all <- NULL
  n_break <- max(length(breakpoint), nbreak)
  for (i in 1:nsim){
    dat_b <- dat[sample.int(n, n, replace = T),]
    res <- fit_pwexp(dat_b$time, dat_b$event, breakpoint, nbreak, max_set, seed=seed+i, trace=FALSE)
    res_all <- rbind(res_all, res)
  }
  if (n_break!=0){
    attr(res_all,'brk') <- res_all[,1:n_break,drop=F]
  }else{
    attr(res_all,'brk') <- NULL
  }
  attr(res_all,'lam') <- res_all[,(n_break+1):(2*n_break+1),drop=F]
  return(res_all)
}

predict_model <- function(event_model_boot, censor_model_boot, n_each=10, train_data, cut, future_enroll, seed=1818, plot=F){
  # future_enroll is a list containing parameters in simdata
  # model is a fitted model
  if (NROW(event_model_boot)!=NROW(censor_model_boot)){
    stop('event_model_boot should have same number of rows as censor_model_boot')
  }
  nsim <- NROW(event_model_boot)
  line_data_fun <- NULL
  line_data_fun_time <- NULL
  ind_cut <- train_data$censor_reason=='cut'
  ind_cut[is.na(ind_cut)] <- FALSE
  for (j in 1:n_each){
    for (i in 1:nsim){
      event_l <- as.numeric(attr(event_model_boot,'lam')[i,])
      censor_l <- as.numeric(attr(censor_model_boot,'lam')[i,])
      event_b <- as.numeric(attr(event_model_boot,'brk')[i,])
      censor_b <- as.numeric(attr(censor_model_boot,'brk')[i,])
      dat_t <- simdata(advanced_dist = list(
        event_dist=function(n)rpwexp(n,event_l,event_b),
        drop_dist=function(n)rpwexp(n,censor_l,censor_b)),
        n_enroll=future_enroll$n_enroll, enroll_rate=future_enroll$enroll_rate, total_sample=future_enroll$total_sample, 
        add_column = c('event','censor_reason','followT_abs','followT'))
      dat_t$followT_abs <- dat_t$followT_abs+cut
      
      dat_t_cut <- simdata(advanced_dist = list(
        event_dist=function(n)rpwexp_conditional(n,train_data$followT[ind_cut], event_l,event_b),
        drop_dist=function(n)rpwexp_conditional(n,train_data$followT[ind_cut], censor_l,censor_b)),
        n_enroll = sum(ind_cut), add_column = c('event','censor_reason','followT_abs','followT'))
      dat_t_cut$followT_abs <- dat_t_cut$followT
      dat_t_cut$followT_abs <- dat_t_cut$followT_abs-train_data$followT[ind_cut]+cut
      dat_pre <- rbind(dat_t, dat_t_cut)
      
      line_data <- plot_event(dat_pre, add=T, additional_event=sum(train_data$event), plot=plot, abs_time = T, col='grey')
      line_data_fun <- c(line_data_fun, suppressWarnings(approxfun(line_data$followT_abs, line_data$n_event, rule=2)))
      line_data_fun_time <- c(line_data_fun_time, suppressWarnings(approxfun(line_data$n_event, line_data$followT_abs, rule=1)))
    }
  }
  return(list(event_fun=line_data_fun,time_fun=line_data_fun_time))
}

predict_rpwexp <- function(predict_model, type='time', at, alpha=0.05){
  #  to predicttype='event' or 'time'
  if (type=='time'){
    pre <- sapply(predict_model$event_fun, function(f)f(at))
  }else {
    pre <- sapply(predict_model$time_fun, function(f)f(at))
  }
  pre[is.na(pre)] <- Inf
  pre <- apply(pre, 1, function(x)quantile(x, c(alpha/2, 0.5, (1-alpha/2))))
  pre <- rbind(at, pre)
  return(t(pre))
}

fit_pwexp_cv <- function(time, event, nfold=5, nsim=100, breakpoint=NULL, nbreak=0, max_set=1000, seed=1818){
  dat <- data.frame(time=time, event=event)
  n <- NROW(dat)
  
  n_break <- max(length(breakpoint), nbreak)
  set.seed(seed)
  cv_like <- NULL
  for (j in 1:nsim){
    ind <- sample(cut(1:n, breaks=nfold, label=FALSE))
    like_inside <- NULL
    for (i in 1:nfold){
      dat_train <- dat[ind!=i,]
      dat_test <- dat[ind==i,]
      md <- fit_pwexp(dat_train$time, dat_train$event, breakpoint, nbreak, max_set, 
                      seed=seed+i+j*nfold, trace=FALSE)
      loglikelihood <- sum(dpwexp(dat_test$time[dat_test$event==1], rate=attr(md,'lam'), breakpoint = attr(md,'brk'), log = T, one_piece = F, safety_check = F))+
        sum(ppwexp(dat_test$time[dat_test$event==0], rate=attr(md,'lam'), lower.tail = F, breakpoint = attr(md,'brk'), log.p = T, one_piece = F, safety_check = F))
      like_inside <- c(like_inside, loglikelihood)
    }
    # like_inside[is.infinite(like_inside)] <- min(like_inside[is.finite(like_inside)])
    cv_like <- c(cv_like, mean(like_inside))
  }
  return(cv_like)
}