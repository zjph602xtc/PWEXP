# simudata function -------------------------------------------------------
# this function generate a simulated clinical trial survival dataset


simdata <- function(group="Group 1", strata='Strata 1', allocation=1, event_lambda=NA, drop_rate=NA,
                    death_lambda=NA, n_rand=NULL, rand_rate=NULL, total_sample=NULL, add_column=c('followT'),
                    simplify=TRUE, advanced_dist=NULL) {
  # total number of subgroups will be '# treatment groups'*'# strata'
  # strata variable will be distributed into each treatment group. For example,
  # group <- c('trt','placebo'), strata <- c('A','B','C'), then there will be six subgroups:
  # trt+A, trt+B, trt+C, placebo+A, placebo+B, placebo+C

  # allocation: the proportion of each subgroup. For example, when group=c('trt','placebo'),
  #             strata=c('A','B','C'), the length of allocation should be 6. Note the value
  #             will be recycled if the length is less than needed.
  # event_lambda: the hazard rate of event. The value will be recycled if the length is less than needed.
  # drop_rate: the drop-out rate (patient/month). The value will be recycled if the
  #            length is less than needed.
  # death_lambda: the hazard rate of death. The value will be recycled if the
  #               length is less than needed.
  # n_rand: (required when rand_rate=NULL) a vector contains the number of randomization each month
  # rand_rate: (required when n_rand=NULL) the randomization rate (patient/month)
  # total_sample: (required when rand_rate is in use) total scheduled sample size
  # add_column: request additional columns
  #           'eventT_abs': absolute event time from the beginning of the trial (=eventT+randT)
  #           'dropT_abs': absolute drop-out time from the beginning of the trial (=dropT+randT)
  #           'deathT_abs': absolute death time from the beginning of the trial (=deathT+randT)
  #           'censor': whether censored
  #           'event': whether having event
  #           'censor_reason': why censored ('drop_out','death','never_event'(eventT=inf))
  #           'followT': follow time (true observed time) from randT
  #           'followT_abs': absolute follow time from the beginning of the trial (=followT+randT)
  # simplify: whether drop unused columns (i.e., the group variable when there is only one group)
  # advanced_dist: use other distributions instead of exponential. The value will be recycled if the
  #               length is less than needed.
  #              e.g., list(event_dist=c(function(n)rexp(n,4), function(n)rweibull(n,3,2), drop_dist=..., death_dist=...)
  #              Each distribution is a function with only one input argument n (sample size)
  #              If any of event_dist, drop_dist, death_dist is missing, then search for event_lambda, drop_rate, death_lambda;
  #              if also missing, then this variable will not be generated

  # output: eventT, dropT, deathT is from randT, not the absolute time from the beginning of the study

  setting <- data.frame(group=rep(group, each=length(strata)), strata=strata, allocation=allocation, event_lambda=event_lambda, drop_rate=drop_rate, death_lambda=death_lambda)
  setting$allocation <- setting$allocation/sum(setting$allocation)
  setting$drop_lambda <- -log(1-setting$drop_rate)   # drop out

  # obtain rand parameter and randT
  if (!is.null(n_rand)){
    Total_d <- sum(n_rand)
    t_rand <- 1:length(n_rand)
    popd <- mapply(runif, n=n_rand, min=c(0,t_rand[-length(t_rand)]), max=t_rand, SIMPLIFY = F)
  }else if (round(rand_rate)==rand_rate){
    n_rand <- rep(rand_rate, ceiling(total_sample/rand_rate))
    n_rand[length(n_rand)] <- total_sample - sum(n_rand[-length(n_rand)])
    Total_d <- sum(n_rand)
    t_rand <- 1:length(n_rand)
    popd <- mapply(runif, n=n_rand, min=c(0,t_rand[-length(t_rand)]), max=t_rand, SIMPLIFY = F)
  }else{
    Total_d <- total_sample
    popd <- runif(total_sample, 0, total_sample/rand_rate)
  }

  # obtain number of subject in each subgroup
  setting$N <- floor(Total_d*setting$allocation)
  ind <- sample.int(NROW(setting), Total_d-sum(setting$N))
  setting$N[ind] <- setting$N[ind]+1
  if (!is.null(advanced_dist)){
    if (any(!names(advanced_dist) %in% c('event_dist','drop_dist','death_dist'))){
      stop('advanced_dist should only contain event_dist, drop_dist, death_dist')
    }
    if (!is.null(advanced_dist$event_dist)){
      if (!inherits(advanced_dist$event_dist,'list')){
        setting$event_lambda <- list(advanced_dist$event_dist)
      }else{
        setting$event_lambda <- advanced_dist$event_dist
      }
    }
    if (!is.null(advanced_dist$drop_dist)){
      if (!inherits(advanced_dist$drop_dist,'list')){
        setting$drop_lambda <- list(advanced_dist$drop_dist)
      }else{
        setting$drop_lambda <- advanced_dist$drop_dist
      }
    }
    if (!is.null(advanced_dist$death_dist)){
      if (!inherits(advanced_dist$death_dist,'list')){
        setting$death_lambda <- list(advanced_dist$death_dist)
      }else{
        setting$death_lambda <- advanced_dist$death_dist
      }
    }
  }

  # generate time
  group_v <- mapply(rep, x=setting$group, each=setting$N, SIMPLIFY = F)
  strata_v <-  mapply(rep, x=setting$strata, each=setting$N, SIMPLIFY = F)
  if (inherits(setting$drop_lambda,c('numeric','logical'))){
    cd <- suppressWarnings(mapply(rexp, n=setting$N, rate=setting$drop_lambda, SIMPLIFY = F))   # censoring time
    try(cd[setting$drop_lambda==0] <- lapply(cd[setting$drop_lambda==0], function(x)rep(Inf, length(x))), silent = T)
  }else{
    cd <- mapply(function(f,n)f(n), f=setting$drop_lambda, n=setting$N, SIMPLIFY = F)
  }
  if (inherits(setting$event_lambda,c('numeric','logical'))){
    ed <- suppressWarnings(mapply(rexp, n=setting$N, rate=setting$event_lambda, SIMPLIFY = F))   # event time
    try(ed[setting$event_lambda==0] <- lapply(ed[setting$event_lambda==0], function(x)rep(Inf, length(x))), silent = T)
  }else{
    ed <- mapply(function(f,n)f(n), f=setting$event_lambda, n=setting$N, SIMPLIFY = F)
  }
  if (inherits(setting$death_lambda,c('numeric','logical'))){
    dd <- suppressWarnings(mapply(rexp, n=setting$N, rate=setting$death_lambda, SIMPLIFY = F))   # death time
    try(dd[setting$death_lambda==0] <- lapply(dd[setting$death_lambda==0], function(x)rep(Inf, length(x))), silent = T)
  }else{
    dd <- mapply(function(f,n)f(n), f=setting$death_lambda, n=setting$N, SIMPLIFY = F)
  }


  DTTE <- data.frame(ID=1:Total_d, group=as.factor(unlist(group_v)), strata=as.factor(unlist(strata_v)), randT=sample(unlist(popd)), eventT=unlist(ed),
                     dropT=unlist(cd), deathT=unlist(dd))


  if (!is.null(add_column)){
    if ('eventT_abs' %in% add_column){
      DTTE$eventT_abs <- DTTE$eventT+DTTE$randT
    }
    if ('deathT_abs' %in% add_column){
      DTTE$deathT_abs <- DTTE$deathT+DTTE$randT
    }
    if ('dropT_abs' %in% add_column){
      DTTE$dropT_abs <- DTTE$dropT+DTTE$randT
    }
    if ('censor_reason' %in% add_column){
      type <- mapply(sum, as.numeric(DTTE$eventT>=DTTE$dropT), 2*(DTTE$eventT>=DTTE$deathT), 4*is.infinite(DTTE$eventT), na.rm = T)
      DTTE$censor_reason <- NA
      DTTE$censor_reason[type==1 | type==5] <- 'drop_out'
      DTTE$censor_reason[type==2 | type==6] <- 'death'
      DTTE$censor_reason[type==3 | type==7] <- ifelse(DTTE$dropT<=DTTE$deathT, 'drop_out','death')[type==3 | type==7]
      DTTE$censor_reason[type==4] <- 'never_event'
    }
    if ('event' %in% add_column | 'censor' %in% add_column){
      event <- mapply(all, DTTE$eventT < DTTE$dropT, DTTE$eventT < DTTE$deathT, DTTE$eventT < Inf, na.rm = T)
      if ('event' %in% add_column){
        DTTE$event <- as.numeric(event)
      }
      if ('censor' %in% add_column){
        DTTE$censor <- as.numeric(!event)
      }
    }
    if ('followT' %in% add_column | 'followT_abs' %in% add_column){
      followT <- mapply(min, DTTE$eventT, DTTE$dropT, DTTE$deathT, na.rm = T)
      if ('followT' %in% add_column){
        DTTE$followT <- followT
      }
      if ('followT_abs' %in% add_column){
        DTTE$followT_abs <- followT+DTTE$randT
      }
    }

  }
  # simplify results
  if (simplify){
    ind <- c('ID', ifelse(length(group)==1,NA,'group'), ifelse(length(strata)==1,NA,'strata'), 'randT', ifelse(all(is.na(setting$event_lambda)), NA, 'eventT'),
             ifelse(all(is.na(setting$drop_lambda)), NA, 'dropT'), ifelse(all(is.na(setting$death_lambda)), NA, 'deathT'), add_column)
    ind <- ind[!is.na(ind)]
    DTTE <- DTTE[,ind]
  }
  rownames(DTTE) <- NULL
  return(DTTE)
}

