
cut_dat <- function(var_enrollT, cut, data, var_followT=NULL, var_followT_abs=NULL, var_censor=NULL,
                    var_event=NULL, var_censor_reason=NULL){
  var <- c(var_followT, var_followT_abs, var_censor, var_event, var_censor_reason)
  if (!is.null(var)){
    notin <- var[!(var %in% colnames(data))]
    if (length(notin)!=0){
      stop(paste0('Variable ',notin, ' is/are not in the column names of the dataset.'))
    }
  }
  # cut is the absolute time from the beginning of the trial
  train <- data[data[var_enrollT] <= cut,,drop=F]
  if (!is.null(var_censor_reason)){
    if (!is.null(var_followT)){
      followT_abs <- train[var_enrollT]+train[var_followT]
      train[var_censor_reason][followT_abs > cut] <- 'cut'
    }else if (!is.null(var_followT_abs)){
      train[var_censor_reason][train[var_followT_abs] > cut] <- 'cut'
    }else{
      stop('var_followT or var_followT_abs must be provided to update var_censor_reason.')
    }
  }
  if (!is.null(var_censor)){
    if (!is.null(var_followT)){
      followT_abs <- train[var_enrollT]+train[var_followT]
      train[var_censor][followT_abs > cut] <- 1
    }else if (!is.null(var_followT_abs)){
      train[var_censor][train[var_followT_abs] > cut] <- 1
    }else{
      stop('var_followT or var_followT_abs must be provided to update var_censor')
    }
  }
  if (!is.null(var_event)){
    if (!is.null(var_followT)){
      followT_abs <- train[var_enrollT]+train[var_followT]
      train[var_event][followT_abs > cut] <- 0
    }else if (!is.null(var_followT_abs)){
      train[var_event][train[var_followT_abs] > cut] <- 0
    }else{
      stop('var_followT or var_followT_abs must be provided to update var_event.')
    }
  }
  if (!is.null(var_followT_abs)){
    train[var_followT_abs][train[var_followT_abs] > cut] <- cut
  }
  if (!is.null(var_followT)){
    followT_abs <- train[var_enrollT]+train[var_followT]
    followT_abs[followT_abs > cut] <- cut
    train[var_followT] <- followT_abs-train[var_enrollT]
  }
  return(train)
}

# tran_censor <- function(dat){
#   dat$event <- ifelse(dat$censor_reason=='censored', 1, 0)
#   dat$event[is.na(dat$event)] <- 0
#   return(dat)
# }
