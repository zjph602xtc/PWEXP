predict.pwexp.fit <- function(event_model, cut_indicator=NULL, analysis_time, censor_model=NULL, type='predictive', n_each=100, future_rand=NULL, seed=1818){
  # future_rand is a list containing parameters in simdata
  # model is a fitted model
  attr(event_model,'lam') <- matrix(attr(event_model,'lam'),nrow=1)
  if (!is.null(attr(event_model,'brk'))){
    attr(event_model,'brk') <- matrix(attr(event_model,'brk'),nrow=1)
  }
  if (!is.null(censor_model)){
    attr(censor_model,'lam') <- matrix(attr(censor_model,'lam'),nrow=1)
    if (!is.null(attr(censor_model,'brk'))){
      attr(censor_model,'brk') <- matrix(attr(censor_model,'brk'),nrow=1)
    }
  }
  res <- predict.boot.pwexp.fit(event_model_boot = event_model, cut_indicator = cut_indicator, analysis_time = analysis_time, censor_model_boot = censor_model, type = type, n_each = n_each, future_rand = future_rand, seed = seed)
  class(res) <- c('predict.pwexp.fit','list')
  return(res)
}

predict.boot.pwexp.fit <- function(event_model_boot, cut_indicator=NULL, analysis_time, censor_model_boot=NULL, type='predictive', n_each=10, future_rand=NULL, seed=1818){
  # future_rand is a list containing parameters in simdata
  # model is a fitted model

  if (!(type %in% c('predictive','confidence'))){
    stop('wrong type argument.')
  }
  if (is.null(censor_model_boot)){
    censormodel <- FALSE
    if (is.null(cut_indicator)){
      warning('Both \'cut_indicator\' and \'censor_model_boot\' are missing. All subjects without events will be regarded as censored at the end of the trial.' )
      cut_indicator <- !attr(event_model_boot, 'para')$event
    }else{
      warning('\'censor_model_boot\' is missing. All subjects without events will be regarded as censored at the end of the trial.')
    }
  }else{
    censormodel <- TRUE
    if (NROW(event_model_boot)!=NROW(censor_model_boot)){
      stop('\'event_model_boot\' should have same number of rows as \'censor_model_boot\'. Please use the same \'nsim\' value when bootstrap both models.')
    }
    if (is.null(cut_indicator)){
      cut_indicator <- !attr(event_model_boot, 'para')$event & !attr(censor_model_boot, 'para')$event
    }
  }
  if (any(is.na(cut_indicator))){
    warning('\'cut_indactor\' contains missing values. All missing values are set to 0')
    cut_indicator[is.na(cut_indicator)] <- 0
  }
  future <- ifelse(is.null(future_rand), FALSE, TRUE)

  censor_l <- 0
  censor_b <- NULL
  nsim <- NROW(event_model_boot)
  line_data_fun <- NULL
  line_data_fun_time <- NULL
  if (type=='confidence'){
    n_each <- 1
  }
  pb <- txtProgressBar(style = 3)
  for (j in 1:n_each){
    if (n_each!=1){
      setTxtProgressBar(pb, j/n_each)
    }
    for (i in 1:nsim){
      event_l <- as.numeric(attr(event_model_boot,'lam')[i,])
      event_b <- as.numeric(attr(event_model_boot,'brk')[i,])
      if (censormodel){
        censor_l <- as.numeric(attr(censor_model_boot,'lam')[i,])
        censor_b <- as.numeric(attr(censor_model_boot,'brk')[i,])
      }
      if (length(event_b)==0) event_b <- NULL
      if (length(censor_b)==0) censor_b <- NULL
      if (future){
        dat_t <- simdata(advanced_dist = list(
          event_dist=function(n)rpwexp(n,event_l,event_b),
          drop_dist=function(n)rpwexp(n,censor_l,censor_b)),
          n_rand=future_rand$n_rand, rand_rate=future_rand$rand_rate, total_sample=future_rand$total_sample,
          add_column = c('event','censor_reason','followT_abs','followT'))
        dat_t$followT_abs <- dat_t$followT_abs+analysis_time
      }else{
        dat_t <- NULL
      }

      para <- attr(event_model_boot,'para')
      dat_t_cut <- simdata(advanced_dist = list(
        event_dist=function(n)rpwexp_conditional(n,para$time[cut_indicator==1], event_l,event_b),
        drop_dist=function(n)rpwexp_conditional(n,para$time[cut_indicator==1], censor_l,censor_b)),
        n_rand = sum(cut_indicator), add_column = c('event','censor_reason','followT_abs','followT'))
      dat_t_cut$followT_abs <- dat_t_cut$followT
      dat_t_cut$followT_abs <- dat_t_cut$followT_abs-para$time[cut_indicator==1]+analysis_time
      dat_pre <- rbind(dat_t, dat_t_cut)

      line_data <- plot_event(time=dat_pre$followT_abs, abs_time = T, add=T, event=dat_pre$event, additional_event=sum(para$event), plot=FALSE, col='grey')
      line_data_fun <- c(line_data_fun, suppressWarnings(approxfun(line_data$time, line_data$n_event, rule=1:2)))
      line_data_fun_time <- c(line_data_fun_time, suppressWarnings(approxfun(line_data$n_event, line_data$time, rule=1)))
    }
  }
  setTxtProgressBar(pb, 1)
  close(pb)
  res <- list(event_fun=line_data_fun,time_fun=line_data_fun_time)
  class(res) <- c('predict.boot.pwexp.fit','list')
  return(res)
}
