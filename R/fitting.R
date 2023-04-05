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

pwexp.fit <- function(time, event, breakpoint=NULL, nbreak=0, max_set=10000, seed=1818, trace=FALSE){
  # n_{fj}/lam_j - sum_{set fj sj}(t_i-d_{j-1}) = n_{set f_j+1 to f_end, s_j+1 to s_end}*(d_j-d_{j-1})
  # f is event data, s is censored data

  # auto search if breakpoint=NULL
  breakpoint <- sort(breakpoint)
  brk_backup <- breakpoint
  time_backup <- time
  event_backup <- event
  ind <- order(time)
  event <- event[ind]
  time <- time[ind]
  time_event <- time[event==1]
  time_noevent <- time[event==0]
  N <- length(time)
  set.seed(seed)
  n_fix_brk <- length(breakpoint)

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

  if (n_fix_brk==0  && nbreak==0){
    lam <- sum(event)/sum(time)
    loglikelihood <- sum(dpwexp(time_event, rate=lam, breakpoint = NULL, log = T, one_piece = T, safety_check = F))+
      sum(ppwexp(time_noevent, rate=lam, lower.tail = F, breakpoint = NULL, log.p = T, one_piece = T, safety_check = F))
    aic <- 2-2*loglikelihood
    bic <- log(length(time))-2*loglikelihood
    res <- data.frame(lam1=lam, likelihood=loglikelihood, AIC=aic, BIC=bic)
    attr(res,'lam') <- lam
    attr(res,'para') <- list(time=time_backup, event=event_backup, breakpoint=brk_backup, nbreak=nbreak)
    class(res) <- c('pwexp.fit', 'data.frame')
    return(res)
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
  attr(res,'para') <- list(time=time_backup, event=event_backup, breakpoint=brk_backup, nbreak=nbreak)
  class(res) <- c('pwexp.fit', 'data.frame')
  return(res)

}

boot.pwexp.fit <- function(x, ...){
  UseMethod("boot.pwexp.fit")
}



boot.pwexp.fit.default <- function(time, event, nsim=100, breakpoint=NULL, nbreak=0, max_set=1000, seed=1818){
  dat <- data.frame(time=time, event=event)
  n <- NROW(dat)
  res_all <- pwexp.fit(dat$time, dat$event, breakpoint, nbreak, max_set, seed=seed, trace=FALSE)

  if (nbreak==0){
    nbreak <- length(attr(res_all, 'lam'))-1
  }
  set.seed(seed)
  pb <- txtProgressBar(style = 3)
  for (i in 1:(nsim-1)){
    setTxtProgressBar(pb, i/(nsim))
    dat_b <- dat[sample.int(n, n, replace = T),]
    res <- pwexp.fit(dat_b$time, dat_b$event, breakpoint, nbreak, max_set, seed=seed+i, trace=FALSE)
    res_all <- rbind(res_all, res)
  }
  setTxtProgressBar(pb, 1)
  close(pb)
  if (nbreak!=0){
    attr(res_all,'brk') <- res_all[,1:nbreak,drop=F]
  }else{
    attr(res_all,'brk') <- NULL
  }
  attr(res_all,'lam') <- res_all[,(nbreak+1):(2*nbreak+1),drop=F]
  class(res_all) <- c('boot.pwexp.fit', 'data.frame')
  return(res_all)
}

boot.pwexp.fit.pwexp.fit <- function(object, nsim=100, max_set=1000, seed=1818){
  para <- attr(object, 'para')
  res <- boot.pwexp.fit.default(time=para$time, event=para$event,
                                nsim=nsim, breakpoint=para$breakpoint, nbreak=para$nbreak,
                                max_set=max_set, seed=seed)
  return(res)
}



cv.pwexp.fit <- function(x, ...){
  UseMethod("cv.pwexp.fit")
}

cv.pwexp.fit.default <- function(time, event, nfold=5, nsim=100, breakpoint=NULL, nbreak=0, max_set=1000, seed=1818){
  dat <- data.frame(time=time, event=event)
  n <- NROW(dat)

  res_all <- pwexp.fit(dat$time, dat$event, breakpoint, nbreak, max_set, seed=seed, trace=FALSE)
  if (nbreak==0){
    nbreak <- length(attr(res_all, 'lam'))-1
  }
  set.seed(seed)
  cv_like <- NULL
  pb <- txtProgressBar(style = 3)
  for (j in 1:nsim){
    setTxtProgressBar(pb, j/nsim)
    ind <- sample(cut(1:n, breaks=nfold, label=FALSE))
    like_inside <- NULL
    for (i in 1:nfold){
      dat_train <- dat[ind!=i,]
      dat_test <- dat[ind==i,]
      md <- pwexp.fit(dat_train$time, dat_train$event, breakpoint, nbreak, max_set,
                      seed=seed+i+j*nfold, trace=FALSE)
      loglikelihood <- sum(dpwexp(dat_test$time[dat_test$event==1], rate=attr(md,'lam'), breakpoint = attr(md,'brk'), log = T, one_piece = F, safety_check = F))+
        sum(ppwexp(dat_test$time[dat_test$event==0], rate=attr(md,'lam'), lower.tail = F, breakpoint = attr(md,'brk'), log.p = T, one_piece = F, safety_check = F))
      like_inside <- c(like_inside, loglikelihood)
    }
    # like_inside[is.infinite(like_inside)] <- min(like_inside[is.finite(like_inside)])
    cv_like <- c(cv_like, mean(like_inside))
  }
  close(pb)
  return(cv_like)
}

cv.pwexp.fit.pwexp.fit <- function(object, nfold=5, nsim=100, max_set=1000, seed=1818){
  para <- attr(object, 'para')
  res <- cv.pwexp.fit.default(time=para$time, event=para$event, nfold=nfold,
                              nsim=nsim, breakpoint=para$breakpoint, nbreak=para$nbreak,
                              max_set=max_set, seed=seed)
  return(res)
}
