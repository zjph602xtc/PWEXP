combine_time <- function(time, tol){
  f_tmp <- time[1]
  f_N <- length(time)
  tol <- (time[f_N] - f_tmp)*tol
  for (i in 2:f_N){
    val <- time[i]
    if ((val - f_tmp)<tol){
      time[i] <- f_tmp
    }else{
      f_tmp <- val
    }
  }
  return(time)
}


get_grid <- function(time, breakpoint, nbreak, max_set=5000, remove_first=TRUE){
  time <- unique(time)
  if (remove_first){
    time <- time[-1]
  }

  if (nbreak==0){
    nbreak <- length(breakpoint)
  }
  if (!is.null(breakpoint)){
    nbreak <- nbreak - length(breakpoint)
    # if (nbreak==0){
    #   return(breakpoint <- matrix(breakpoint, nrow=1))
    # }
    # if (nbreak<0){
    #   stop('nbreak is the total number of breakpoints, which must be equal or larger than the length of \'breakpoint\'')
    # }
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
  if (length(time) < nbreak){
    setting <- matrix(-Inf, ncol=nbreak, nrow=1)
  }else{
    setting <- t(combn(time, nbreak))
  }


  if (NROW(setting) > max_set){
    setting <- setting[sort(sample.int(NROW(setting), max_set, replace = F)),,drop=F]
  }

  if (!is.null(breakpoint)){
    setting <- t(apply(setting, 1, function(x)sort(c(breakpoint, x))))
  }
  return(setting)
}



pwexp.fit <- function(time, event, breakpoint=NULL, nbreak=0, exclude_int=NULL, min_pt_tail=5, max_set=10000, seed=1818, trace=FALSE, optimizer='mle', tol=1e-4){
  # n_{fj}/lam_j - sum_{set fj sj}(t_i-d_{j-1}) = n_{set f_j+1 to f_end, s_j+1 to s_end}*(d_j-d_{j-1})
  # f is event data, s is censored data

  # auto search if breakpoint=NULL
  if (!(optimizer %in% c('hybrid', 'ols', 'mle'))){
    stop('wrong optimizer!')
  }
  if (optimizer != 'mle' & !is.null(exclude_int)){
    warning('\'exclude_int\' only works for optimizer=\'mle\'.')
  }
  if (optimizer != 'mle' & min_pt_tail!=5){
    warning('\'min_pt_tail\' only works for optimizer=\'mle\'.')
  }

  event <- as.numeric(event)
  breakpoint <- sort(breakpoint)
  brk_backup <- breakpoint
  time_backup <- time
  event_backup <- event
  ind <- order(time)
  event <- event[ind]
  time <- time[ind]
  if (tol!=0){
    time <- combine_time(time, tol)
  }
  time_event <- time[event==1]
  time_noevent <- time[event==0]
  N <- length(time)
  set.seed(seed)
  n_fix_brk <- length(breakpoint)

  if (n_fix_brk==0  && is.null(nbreak)){
    nbreak <- ceiling(8*length(time_event)^0.2)
    message(paste0('Number of breakpoints = ', nbreak))
  }

  if (n_fix_brk==0  && nbreak==0){
    lam <- sum(event)/sum(time)
    loglikelihood <- sum(PWEXP::dpwexp(time_event, rate=lam, breakpoint = NULL, log = T, one_piece = T, safety_check = F))+
      sum(PWEXP::ppwexp(time_noevent, rate=lam, lower.tail = F, breakpoint = NULL, log.p = T, one_piece = T, safety_check = F))
    aic <- 2-2*loglikelihood
    bic <- log(length(time))-2*loglikelihood
    if (lam==0){
      lam <- loglikelihood <- -Inf
      aic <- bic <- Inf
      warning('Incorrect result returned. Please check the total number of events is not 0')
    }
    res <- data.frame(lam1=lam, likelihood=loglikelihood, AIC=aic, BIC=bic)
    attr(res,'lam') <- lam
    attr(res,'para') <- list(time=time_backup, event=event_backup, breakpoint=brk_backup, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail)
    class(res) <- c('pwexp.fit', 'data.frame')
    return(res)
  }

  while (n_fix_brk > 0){
    tmpi <- findInterval(time, vec=c(0, breakpoint, Inf))
    numerator <- ctapply(event, tmpi, sum)
    if (all(numerator > 0) && all((1:(n_fix_brk+1)) %in% names(numerator))){
      break
    }
    ind <- min(which(numerator==0), setdiff(1:(n_fix_brk+1), names(numerator)), na.rm=T)
    if (is.na(breakpoint[ind])){
      warning(paste0(breakpoint[ind-1],' is too large. No event after this breakpoint. We will remove it.'))
      breakpoint <- breakpoint[-(ind-1)]
    }else if (ind==1){
      warning(paste0(breakpoint[ind],' is too early. No event before this breakpoint. We will remove it.'))
      breakpoint <- breakpoint[-ind]
    }else{
      warning(paste0('There are NO events between ', breakpoint[ind-1], ' and ',
                     breakpoint[ind], '. We will use their mid point ',
                     mean(breakpoint[(ind-1):ind]), ' as a breakpoint. '))
      breakpoint[ind-1] <- mean(breakpoint[(ind-1):ind])
      breakpoint <- breakpoint[-ind]
    }
    n_fix_brk <- length(breakpoint)
    if (n_fix_brk==0){
      breakpoint <- NULL
    }
  }


  if (!is.null(breakpoint)){
    if (nbreak==0 | nbreak == length(breakpoint)){
      breakpoint <- matrix(breakpoint, nrow=1)
      optimizer <- 'fixed'
    }
    if (nbreak!=0 & nbreak < length(breakpoint)){
      stop('nbreak is the total number of breakpoints, which must be equal or larger than the length of \'breakpoint\'')
    }
  }

  if (optimizer=='hybrid' | optimizer=='ols' ){
    Sfun <- survfit(Surv(time, event)~1)
    Sfun <- data.frame(x=Sfun$time, y=log(Sfun$surv))
    Sfun <- na.omit(Sfun)
    Sfun <- Sfun[!is.infinite(Sfun$y),]
    invisible(capture.output(seg_brk <- try(segmented::segmented(lm(y~x, data=Sfun), npsi = nbreak-n_fix_brk, fixed.psi = breakpoint)$psi, silent=TRUE)))
    if (is.null(seg_brk) | any(class(seg_brk)=='try-error')){
      optimizer <- 'mle'
    }

    if (optimizer=='ols'){
      add_by_ols <- c(rep(TRUE, nbreak-n_fix_brk), rep(FALSE, n_fix_brk))
      tmp <- order(c(seg_brk[,2], breakpoint))
      add_by_ols <- add_by_ols[tmp]
      breakpoint <- matrix(c(seg_brk[,2], breakpoint)[tmp], nrow=1)

      tmp_n_fix_brk <- length(breakpoint)
      while (tmp_n_fix_brk > 0){
        tmpi <- findInterval(time, vec=c(0, breakpoint, Inf))
        numerator <- ctapply(event, tmpi, sum)
        if (all(numerator > 0) && all((1:(tmp_n_fix_brk+1)) %in% names(numerator))){
          break
        }
        else{
          optimizer <- 'mle'
        }
        ind <- min(which(numerator==0), setdiff(1:(tmp_n_fix_brk+1), names(numerator)), na.rm=T)
        if (is.na(breakpoint[ind])){
          breakpoint <- breakpoint[-(ind-1)]
          add_by_ols <- add_by_ols[-(ind-1)]
        }else if (ind==1){
          breakpoint <- breakpoint[-ind]
          add_by_ols <- add_by_ols[-ind]
        }else{
          breakpoint[ind-1] <- mean(breakpoint[(ind-1):ind][add_by_ols[(ind-1):ind]])
          breakpoint <- breakpoint[-ind]
          add_by_ols <- add_by_ols[-ind]
        }
        tmp_n_fix_brk <- length(breakpoint)
        if (tmp_n_fix_brk==0){
          breakpoint <- NULL
        }
      }
    }else if (optimizer=='hybrid'){
      candidate_time <- unlist(mapply(function(lb,ub)time[time >=lb & time<=ub], seg_brk[,2]-2*seg_brk[,3], seg_brk[,2]+2*seg_brk[,3]))
      if (length(candidate_time) <= 2*(nbreak-length(breakpoint))){
        optimizer <- 'mle'
      }else{
        breakpoint <- get_grid(sort(candidate_time), breakpoint, nbreak, max_set, remove_first=FALSE)
      }
    }
  }

  if (optimizer=='mle'){
    change_cand <- FALSE
    if(min_pt_tail!=0 | min_pt_tail!=1){
      change_cand <- TRUE
      if ((length(time_event)-min_pt_tail) < 2){
        warning('Insufficient data due to too many points reserved for tail. Ignore \'min_pt_tail\'.')
        candidate_time <- time
      }else{
        candidate_time <- time[time <= time_event[length(time_event)-min_pt_tail]]
        if (length(candidate_time) <= 2*(nbreak-length(breakpoint))){
          warning('Insufficient data due to too many points reserved for tail. Ignore \'min_pt_tail\'.')
          candidate_time <- time
        }
      }
    }
    if(!is.null(exclude_int)){
      change_cand <- TRUE
      candidate_time <- time[time < exclude_int[1] | time > exclude_int[2]]
      if (length(candidate_time) <= 2*(nbreak-length(breakpoint))){
        warning('Insufficient data due to wide coverage of \'exclude_int\'. Ignore \'exclude_int\'.')
        candidate_time <- time
      }
    }
    if (!change_cand){
      candidate_time <- time
    }
    breakpoint <- get_grid(candidate_time, breakpoint, nbreak, max_set)
  }


  if (is.infinite(breakpoint[1,1])){
    res <- matrix(-Inf, ncol = 2*NCOL(breakpoint)+2, nrow=1)
  }else{
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
      loglikelihood <- sum(PWEXP::dpwexp(time_event, rate=lam, breakpoint = brk0, log = T, one_piece = F, safety_check = F))+
        sum(PWEXP::ppwexp(time_noevent, rate=lam, lower.tail = F, breakpoint = brk0, log.p = T, one_piece = F, safety_check = F))
      res[i,] <- c(brk0, lam, loglikelihood)
    }
  }

  res <- data.frame(res)
  if (!trace){
    res <- res[which.max(res[,NCOL(res)])[1],,drop=F]
  }
  n_k <- 2*max(n_fix_brk, nbreak)+1
  aic <- 2*(n_k-n_fix_brk)-2*res[,NCOL(res)]
  bic <- (n_k-n_fix_brk)*log(N)-2*res[,NCOL(res)]
  res <- cbind(res, aic, bic)
  if (!trace){
    attr(res,'lam') <- as.numeric(res[,(NCOL(breakpoint)+1):(2*NCOL(breakpoint)+1)])
    attr(res,'brk') <- as.numeric(res[,1:NCOL(breakpoint)])
  }
  colnames(res) <- c(paste0('brk', 1:NCOL(breakpoint)), paste0('lam', 1:(NCOL(breakpoint)+1)), 'likelihood','AIC','BIC')
  attr(res,'para') <- list(time=time_backup, event=event_backup, breakpoint=brk_backup, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail)
  class(res) <- c('pwexp.fit', 'data.frame')
  if (any(is.infinite(as.numeric(res)))){
    warning('Incorrect result returned. Please check the number of events is at least 2 more than the number of breakpoints. ')
  }
  return(res)
}

boot.pwexp.fit <- function(time, ...){
  UseMethod("boot.pwexp.fit")
}



boot.pwexp.fit.default <- function(time, event, nsim=100, breakpoint=NULL, nbreak=0, exclude_int=NULL, min_pt_tail=5, max_set=1000, seed=1818, optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...){
  dat <- data.frame(time=time, event=event)
  n <- NROW(dat)
  res_all <- pwexp.fit(time=dat$time, event=dat$event, breakpoint=breakpoint, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, seed=seed, trace=FALSE, optimizer=optimizer, tol=tol)

  ind <- order(dat$time)
  dat <- dat[ind,,drop=F]
  if (tol!=0){
    dat$time <- combine_time(dat$time, tol)
  }

  if (nbreak==0){
    nbreak <- length(attr(res_all, 'lam'))-1
  }
  pb <- txtProgressBar(max = nsim, style = 3)

  if (parallel){
    doSNOW::registerDoSNOW(cl <- parallel::makeCluster(mc.core))
    `%dopar%` <- foreach::`%dopar%`
    res_all_tp <- foreach::foreach(i=1:(nsim-1), .combine = 'rbind', .inorder = FALSE, .errorhandling = 'remove', .packages = 'PWEXP', .options.snow=list(progress=function(n)setTxtProgressBar(pb, n))) %dopar% {
      dat_b <- dat[sample.int(n, n, replace = T),]
      res <- suppressWarnings(pwexp.fit(time=dat_b$time, event=dat_b$event, breakpoint=breakpoint, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, seed=seed+i, trace=FALSE, optimizer=optimizer, tol=0))
    }
    res_all <- rbind(res_all, res_all_tp)
    parallel::stopCluster(cl)
  }else{
    for (i in 1:(nsim-1)){
      setTxtProgressBar(pb, i)
      dat_b <- dat[sample.int(n, n, replace = T),]
      res <- suppressWarnings(pwexp.fit(time=dat_b$time, event=dat_b$event, breakpoint=breakpoint, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, seed=seed+i, trace=FALSE, optimizer=optimizer, tol=0))
      res_all <- rbind(res_all, res)
    }
  }

  setTxtProgressBar(pb, nsim)
  close(pb)
  res_all[is.infinite(res_all[,1]),] <- NA
  res_all[is.na(res_all[,1]),] <- suppressWarnings(matrix(colMeans(res_all, na.rm=T), ncol=NCOL(res_all), nrow=sum(is.na(res_all[,1])), byrow = T))
  if (nbreak!=0){
    attr(res_all,'brk') <- res_all[,1:nbreak,drop=F]
  }else{
    attr(res_all,'brk') <- NULL
  }
  attr(res_all,'lam') <- res_all[,(nbreak+1):(2*nbreak+1),drop=F]
  class(res_all) <- c('boot.pwexp.fit', 'data.frame')
  return(res_all)
}

boot.pwexp.fit.pwexp.fit <- function(time, nsim=100, max_set=1000, seed=1818, optimizer='mle', tol=1e-4,
                                     parallel=FALSE, mc.core=4, ...){
  object <- time
  para <- attr(object, 'para')
  res <- boot.pwexp.fit.default(time=para$time, event=para$event,
                                nsim=max(1,nsim-1), breakpoint=para$breakpoint, nbreak=para$nbreak,
                                exclude_int=para$exclude_int, min_pt_tail=para$min_pt_tail,
                                max_set=max_set, seed=seed, optimizer=optimizer, tol=tol,
                                parallel=parallel, mc.core=mc.core)
  res_combined <- rbind(object, res)
  attr(res_combined, 'brk') <- rbind(attr(object, 'brk'), attr(res, 'brk'))
  attr(res_combined, 'lam') <- rbind(attr(object, 'lam'), attr(res, 'lam'))
  class(res_combined) <- c('boot.pwexp.fit', 'data.frame')
  return(res_combined)
}



cv.pwexp.fit <- function(time, ...){
  UseMethod("cv.pwexp.fit")
}

cv.pwexp.fit.default <- function(time, event, nfold=5, nsim=100, breakpoint=NULL, nbreak=0, exclude_int=NULL, min_pt_tail=5, max_set=1000, seed=1818, optimizer='mle', tol=1e-4, parallel=FALSE, mc.core=4, ...){
  dat <- data.frame(time=time, event=event)
  n <- NROW(dat)
  res_all <- pwexp.fit(time=dat$time, event=dat$event, breakpoint=breakpoint, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, seed=seed, trace=FALSE, optimizer=optimizer, tol=tol)
  if (nbreak==0){
    nbreak <- length(attr(res_all, 'lam'))-1
  }

  ind <- order(dat$time)
  dat <- dat[ind,,drop=F]
  if (tol!=0){
    dat$time <- combine_time(dat$time, tol)
  }

  cv_like <- NULL
  pb <- txtProgressBar(max=nsim, style = 3)

  if (parallel){
    doSNOW::registerDoSNOW(cl <- parallel::makeCluster(mc.core))
    `%dopar%` <- foreach::`%dopar%`
    cv_like <- foreach::foreach(j=1:nsim, .combine = 'c', .inorder = FALSE, .errorhandling = 'remove', .packages = 'PWEXP', .options.snow=list(progress=function(n)setTxtProgressBar(pb, n))) %dopar% {
      ind <- sample(cut(1:n, breaks=nfold, label=FALSE))
      like_inside <- NULL
      for (i in 1:nfold){
        dat_train <- dat[ind!=i,]
        dat_test <- dat[ind==i,]
        md <- pwexp.fit(time=dat_train$time, event=dat_train$event, breakpoint=breakpoint, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, seed=seed+i+j*nfold, trace=FALSE, optimizer=optimizer, tol=0)
        if (is.infinite(md[1,1])){
          next
        }
        loglikelihood <- sum(PWEXP::dpwexp(dat_test$time[dat_test$event==1], rate=attr(md,'lam'), breakpoint = attr(md,'brk'), log = T, one_piece = F, safety_check = F))+
          sum(PWEXP::ppwexp(dat_test$time[dat_test$event==0], rate=attr(md,'lam'), lower.tail = F, breakpoint = attr(md,'brk'), log.p = T, one_piece = F, safety_check = F))
        like_inside <- c(like_inside, loglikelihood)
      }
      # like_inside[is.infinite(like_inside)] <- min(like_inside[is.finite(like_inside)])
      mean(like_inside)
    }
    parallel::stopCluster(cl)
  }else{
    for (j in 1:nsim){
      setTxtProgressBar(pb, j)
      ind <- sample(cut(1:n, breaks=nfold, label=FALSE))
      like_inside <- NULL
      for (i in 1:nfold){
        dat_train <- dat[ind!=i,]
        dat_test <- dat[ind==i,]
        md <- pwexp.fit(time=dat_train$time, event=dat_train$event, breakpoint=breakpoint, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, seed=seed+i+j*nfold, trace=FALSE, optimizer=optimizer, tol=0)
        if (is.infinite(md[1,1])){
          next
        }
        loglikelihood <- sum(PWEXP::dpwexp(dat_test$time[dat_test$event==1], rate=attr(md,'lam'), breakpoint = attr(md,'brk'), log = T, one_piece = F, safety_check = F))+
          sum(PWEXP::ppwexp(dat_test$time[dat_test$event==0], rate=attr(md,'lam'), lower.tail = F, breakpoint = attr(md,'brk'), log.p = T, one_piece = F, safety_check = F))
        like_inside <- c(like_inside, loglikelihood)
      }
      # like_inside[is.infinite(like_inside)] <- min(like_inside[is.finite(like_inside)])
      cv_like <- c(cv_like, mean(like_inside))
    }
  }

  close(pb)
  cv_like <- cv_like[!is.na(cv_like)]
  return(cv_like)
}

cv.pwexp.fit.pwexp.fit <- function(time, nfold=5, nsim=100, max_set=1000, seed=1818, optimizer='mle', tol=1e-4,
                                   parallel=FALSE, mc.core=4, ...){
  object <- time
  para <- attr(object, 'para')
  res <- cv.pwexp.fit.default(time=para$time, event=para$event, nfold=nfold,
                              nsim=nsim, breakpoint=para$breakpoint, nbreak=para$nbreak, exclude_int=para$exclude_int, min_pt_tail=para$min_pt_tail,
                              max_set=max_set, seed=seed, optimizer=optimizer, tol=tol, parallel=parallel, mc.core=mc.core)
  return(res)
}
