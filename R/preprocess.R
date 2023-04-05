
cut_dat <- function(cut, data, var_enrollT=NULL, var_followT=NULL, var_followT_abs=NULL, var_censor=NULL,
                    var_event=NULL, var_censor_reason='status_at_end'){
  var <- c(var_followT, var_followT_abs, var_censor, var_event)
  if (!is.null(var)){
    notin <- var[!(var %in% colnames(data))]
    if (length(notin)!=0){
      stop(paste0('Variable ',notin, ' is/are not in the column names of the dataset.'))
    }
  }
  if (!(var_censor_reason %in% colnames(data))){
    warning(paste0('Censor reason variable ',var_censor_reason, ' is not in the column names of the dataset. We will create this variable.'))
    data[,var_censor_reason] <- NA
  }
  # cut is the absolute time from the beginning of the trial
  var_enrollT_missing <- FALSE
  if (is.null(var_enrollT)){
    var_enrollT <- paste0(sample(LETTERS, 64, replace = T), collapse = '')
    data[,var_enrollT] <- 0
    var_enrollT_missing <- TRUE
    warning('\'var_enrollT\' is missing. Treat all subjects enrolled at time=0.')
  }
  train <- data[data[var_enrollT] <= cut,,drop=F]
  if (!is.null(var_censor_reason)){
    if (!is.null(var_followT)){
      followT_abs <- train[var_enrollT]+train[var_followT]
      train[var_censor_reason][followT_abs > cut] <- 'cut'
    }else if (!is.null(var_followT_abs)){
      train[var_censor_reason][train[var_followT_abs] > cut] <- 'cut'
    }else{
      stop('\'var_followT\' or \'var_followT_abs\' must be provided to update var_censor_reason.')
    }
  }
  if (!is.null(var_censor)){
    if (!is.null(var_followT)){
      followT_abs <- train[var_enrollT]+train[var_followT]
      train[var_censor][followT_abs > cut] <- 1
    }else if (!is.null(var_followT_abs)){
      train[var_censor][train[var_followT_abs] > cut] <- 1
    }else{
      stop('\'var_followT\' or \'var_followT_abs\' must be provided to update var_censor')
    }
  }
  if (!is.null(var_event)){
    if (!is.null(var_followT)){
      followT_abs <- train[var_enrollT]+train[var_followT]
      train[var_event][followT_abs > cut] <- 0
    }else if (!is.null(var_followT_abs)){
      train[var_event][train[var_followT_abs] > cut] <- 0
    }else{
      stop('\'var_followT\' or \'var_followT_abs\' must be provided to update var_event.')
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
  if (var_enrollT_missing){
    train <- train[, !(colnames(train)==var_enrollT)]
  }
  return(train)
}

# tran_censor <- function(dat){
#   dat$event <- ifelse(dat$censor_reason=='censored', 1, 0)
#   dat$event[is.na(dat$event)] <- 0
#   return(dat)
# }
