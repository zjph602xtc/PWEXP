sim_followup <- function(at, type = 'calander', group="Group 1", strata='Strata 1', allocation=1,
                         event_lambda=NA, drop_rate=NA, death_lambda=NA, n_rand=NULL, rand_rate=NULL,
                         total_sample=NULL, min_follow=0, by_group=FALSE, by_strata=FALSE,
                         advanced_dist=NULL, stat=c(mean, median, sum),
                         count_in_min_follow=FALSE, count_insufficient_event=FALSE, start_date=NULL, rep=300, seed=1818){
  # this function simulates the average followup time at IAs

  # at: calculate the average follow-up time at 'at' randomized subjects (type='sample') or at 'at' events (type='event')
  #     or at time 'at' (type='calender')
  # type: 'sample'/'event'/'calendar'.  calculate at the number of randomized subjects or events, or at a fixed calendar time
  # rand_rate: (required when n_rand=NULL) the randomization rate (patient/month)
  # n_rand: (required when rand_rate=NULL) a vector contains the number of randomization each month
  # min_follow: (optional) minimum follow up time for the last subject at time 'at'
  # total_sample: (required when type='event' & n_rand=NULL) total scheduled sample size

  # group, strata, allocation: see help of simdata function
  # event_lambda: (required when type='event') the hazard rate of event
  # drop_rate: (optional) the drop-out rate (patient/month)
  # death_lambda: (optional) the hazard rate of death
  # advanced_dist: (optional) see help of simdata function

  # by_group: whether show results stratified by group variable
  # by_strata: whether show results stratified by strata variable
  # stat: which stat is calculated for the follow up time. Can be a user defined function that takes a vector as input and
  #       returns a single value

  # count_in_min_follow: whether count subjects who are randomized after (time of 'at') but before (time of 'at' + min_follow)
  # count_insufficient_event: the method to deal with the case where total event number never reaches required number of events in 'at'.
  #                           If FALSE, then skip the case and show a warning; If TRUE, then use the time of end of the study (the time when
  #                           all subjects die or are censored or have events).
  # start_data: the start date of the trial, in the format: "2000-01-30"
  # rep: number of simulations
  # seed: the random seed



  # function begins ---------------------------------------------------------
  if (!is.null(rand_rate) & type=='sample'){
    total_sample <- max(at)
  }
  if (!(type %in% c('sample','event','calendar'))){
    stop('wrong \'type\'.')
  }
  stat_name <- as.character(substitute(stat))
  stat_name <- setdiff(stat_name,'c')
  has_event <- all(!is.na(event_lambda)) | !is.null(advanced_dist$event_dist)

  set.seed(seed)
  T_all <- NULL
  T_by_group <- tp_by_group <- NULL
  T_by_strata <- tp_by_strata <- NULL
  T_by_group_strata <- tp_by_group_strata <- NULL
  for (iter in 1:rep){
    # simulate data set
    dat <- simdata(group, strata, allocation, event_lambda, drop_rate,
                   death_lambda, n_rand, rand_rate, total_sample, NULL, simplify=F, advanced_dist)
    # when event parameter exists
    if (has_event){
      dat$event <- mapply(all, dat$eventT < dat$dropT, dat$eventT < dat$deathT, dat$eventT < Inf, na.rm = T)
      dat$eventT_abs <- dat$randT+dat$eventT
      dat <- dat[order(dat$eventT_abs),]
      dat$cumevent <- cumsum(dat$event)
    }
    for (i in 1:length(at)){
      # calculate cut-off time (analysis time)
      if (type=='sample'){
        # type: sample size
        dat <- dat[order(dat$randT),]
        Cut_T1 <- dat$randT[at[i]] + min_follow
      }else if (type=='calendar'){
        # type: calendar time
        Cut_T1 <- at[i] + min_follow
      }else if (type=='event'){
        # type: event number
        Cut_T1 <- dat$eventT_abs[dat$cumevent==at[i]][1] + min_follow
        if (is.na(Cut_T1) & count_insufficient_event){
          Cut_T1 <- max(dat$randT+ mapply(min, dat$eventT, dat$dropT, dat$deathT, na.rm=T), na.rm = T) + min_follow
        }else if (is.na(Cut_T1) & !count_insufficient_event){
          warning(paste0('In iteration ', iter, ', the total number of events does not reach ',at[i], ', so we skip this iteration.'))
          next
        }
      }
      # get subset of subjects used to calculate follow up time
      if (count_in_min_follow){
        tmp <- dat[dat$randT <= Cut_T1, ]
      }else{
        tmp <- dat[dat$randT <= (Cut_T1 - min_follow), ]
      }
      # calculate follow-up time
      tmp$followT <- mapply(min, Cut_T1-tmp$randT, tmp$dropT, tmp$deathT, na.rm=T)

      # for the whole dataset without stratification
      tp_all <- aggregate(followT~1, data=tmp, FUN=function(x)sapply(c(length,stat), function(f)f(x)[1]),drop=F)
      tp_all <- cbind(at=at[i], analysis_time=Cut_T1, tp_all)
      if (has_event){
        tmpevent <- tmp[tmp$eventT_abs<=Cut_T1,]
        tp_all <- cbind(tp_all, aggregate(event~1, data=tmpevent, FUN = sum, drop=F))
      }
      T_all <- rbind(T_all, tp_all)
      if (by_group){
        tp_by_group <- aggregate(followT~group, data=tmp, FUN=function(x)sapply(c(length,stat), function(f)f(x)[1]),drop=F)
        tp_by_group <- cbind(at=at[i], analysis_time=Cut_T1, tp_by_group)
        if (has_event){
          tp_by_group <- cbind(tp_by_group, aggregate(event~group, data=tmpevent, FUN = sum, drop=F))
        }
        T_by_group <- rbind(T_by_group, tp_by_group)
      }
      if (by_strata){
        tp_by_strata <- aggregate(followT~strata, data=tmp, FUN=function(x)sapply(c(length,stat), function(f)f(x)[1]),drop=F)
        tp_by_strata <- cbind(at=at[i], analysis_time=Cut_T1, tp_by_strata)
        if (has_event){
          tp_by_strata <- cbind(tp_by_strata, aggregate(event~strata, data=tmpevent, FUN = sum, drop=F))
        }
        T_by_strata <- rbind(T_by_strata, tp_by_strata)
      }
      if (all(by_group, by_strata)){
        tp_by_group_strata <- aggregate(followT~group+strata, data=tmp, FUN=function(x)sapply(c(length,stat), function(f)f(x)[1]),drop=F)
        tp_by_group_strata <- cbind(at=at[i], analysis_time=Cut_T1, tp_by_group_strata)
        if (has_event){
          tp_by_group_strata <- cbind(tp_by_group_strata, aggregate(event~group+strata, data=tmpevent, FUN = sum, drop=F))
        }
        T_by_group_strata <- rbind(T_by_group_strata, tp_by_group_strata)
      }
    }
  }

  all_res_var <- c('at','group','strata','analysis_time', 'analysis_time_c', 'event')
  cluster <- c('group','strata')
  all_res <- list(T_all=T_all, T_by_group=T_by_group, T_by_strata=T_by_strata, T_by_group_strata=T_by_group_strata)
  all_res <- all_res[!sapply(all_res,is.null)]
  for (i in 1:length(all_res)){
    tmp <- all_res[[i]]
    tmp <- cbind(tmp[,all_res_var[all_res_var %in% colnames(tmp)]], tmp$followT)
    colnames(tmp) <-  c(all_res_var[all_res_var %in% colnames(tmp)],'subjects',stat_name)
    if (!is.null(start_date)) tmp$analysis_time_c <- as.Date(start_date)+tmp$analysis_time*30.4375
    tmp <- suppressWarnings(aggregate(tmp, as.list(tmp[,c('at',cluster[cluster %in% colnames(tmp)]),drop=F]), mean))
    tmp <- tmp[order(tmp$at),]
    all_res[[i]] <- tmp[,c(all_res_var[all_res_var %in% colnames(tmp)],'subjects',stat_name)]
  }
  return(all_res)
}


