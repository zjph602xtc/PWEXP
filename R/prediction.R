predict.pwexpm <- function(object, cut_indicator=NULL, analysis_time, censor_model=NULL, n_each=100, future_rand=NULL, seed=1818, ...){
  # future_rand is a list containing parameters in simdata
  # model is a fitted model
  event_model <- object
  res <- predict.boot.pwexpm(object = event_model, cut_indicator = cut_indicator, analysis_time = analysis_time, censor_model = censor_model, n_each = n_each, future_rand = future_rand, seed = seed)
  class(res) <- c('predict.pwexpm','list')
  return(res)
}

predict.boot.pwexpm <- function(object, cut_indicator=NULL, analysis_time, censor_model=NULL, n_each=10, future_rand=NULL, seed=1818, ...){
  # future_rand is a list containing parameters in simdata
  # model is a fitted model
  event_model_boot <- object
  if (!is.null(seed)) set.seed(seed)
  if (is.null(censor_model)){
    censormodel <- FALSE
    if (is.null(cut_indicator)){
      warning('Both \'cut_indicator\' and \'censor_model\' are missing. All subjects without events will be regarded as censored at the end of the trial.' )
      cut_indicator <- !event_model_boot$para$event
    }else{
      warning('\'censor_model\' is missing. All subjects without events will be regarded as censored at the end of the trial.')
    }
  }else{
    censormodel <- TRUE
    if (NROW(event_model_boot$lam)!=NROW(censor_model$lam)){
      stop('\'event_model_boot\' should have same number of rows as \'censor_model\'. Please use the same \'nsim\' value when bootstrap both models.')
    }
    if (is.null(cut_indicator)){
      cut_indicator <- !event_model_boot$para$event & !censor_model$para$event
    }
  }
  if (any(is.na(cut_indicator))){
    warning('\'cut_indactor\' contains missing values. All missing values are set to 0')
    cut_indicator[is.na(cut_indicator)] <- 0
  }
  future <- ifelse(is.null(future_rand), FALSE, TRUE)
  if (future){
    if (!is.null(future_rand$total_sample)){
      if (future_rand$total_sample <=0){
        future <- FALSE
      }
    }
    if (!is.null(future_rand$n_rand)){
      if (sum(future_rand$n_rand) <= 0){
        future <- FALSE
      }
    }
  }

  censor_l <- 0
  censor_b <- NULL
  nsim <- NROW(event_model_boot$lam)
  line_data_fun <- list()
  # if (type=='confidence'){
  #   n_each <- 1
  # }
  pb <- txtProgressBar(style = 3)
  for (j in 1:n_each){
    if (n_each!=1){
      setTxtProgressBar(pb, j/n_each)
    }
    for (i in 1:nsim){
      event_l <- unlist(event_model_boot$lam[i,])
      event_b <- unlist(event_model_boot$brk[i,])
      if (censormodel){
        censor_l <- unlist(censor_model$lam[i,])
        censor_b <- unlist(censor_model$brk[i,])
      }
      if (length(event_b)==0) event_b <- NULL
      if (length(censor_b)==0) censor_b <- NULL
      if (future){
        dat_t <- simdata(advanced_dist = list(
          event_dist=function(n)PwePred::rpwexpm(n,event_l,event_b),
          drop_dist=function(n)PwePred::rpwexpm(n,censor_l,censor_b)),
          n_rand=future_rand$n_rand, rand_rate=future_rand$rand_rate, total_sample=future_rand$total_sample,
          add_column = c('event','censor_reason','followT_abs','followT'))
        dat_t$followT_abs <- dat_t$followT_abs+analysis_time
      }else{
        dat_t <- NULL
      }

      para <- event_model_boot$para
      dat_t_cut <- simdata(advanced_dist = list(
        event_dist=function(n)rpwexpm_conditional(n,para$time[cut_indicator==1], event_l,event_b),
        drop_dist=function(n)rpwexpm_conditional(n,para$time[cut_indicator==1], censor_l,censor_b)),
        n_rand = sum(cut_indicator), add_column = c('event','censor_reason','followT_abs','followT'))
      dat_t_cut$followT_abs <- dat_t_cut$followT
      dat_t_cut$followT_abs <- dat_t_cut$followT_abs-para$time[cut_indicator==1]+analysis_time
      dat_pre <- rbind(dat_t, dat_t_cut)

      line_data <- plot_event(time=dat_pre$followT_abs, abs_time = T, add=T, event=dat_pre$event, additional_event=sum(para$event), plot=FALSE, col='grey')
      line_data <- rbind(c(analysis_time, sum(para$event)), line_data)
      flag <- tryCatch({
        tmp_line_dat_fun <- suppressWarnings(approxfun(line_data$time, line_data$n_event, rule=1:2))
      }, error = function(e){
        e
      })
      if (inherits(flag, "error")){
        next
      }else{
        line_data_fun[[paste0('boot_sample_',i , '_n_each_', j)]] <- tmp_line_dat_fun
        # line_data_fun <- c(line_data_fun,tmp_line_dat_fun)
      }
    }
  }
  setTxtProgressBar(pb, 1)
  close(pb)
  res <- list(event_fun=line_data_fun)
  res$event_model <- object
  res$censor_medol <- censor_model
  res$nsim <- nsim
  res$bootstrap <- ifelse(nsim > 1, TRUE, FALSE)
  res$para <- list(cut_indicator=cut_indicator, analysis_time=analysis_time, n_each=n_each, future_rand=future_rand, seed=seed)
  class(res) <- c('predict.boot.pwexpm','list')
  return(res)
}
