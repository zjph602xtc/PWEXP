coef.pwexpm <- function(object, parm = 'all', ...){
  if (parm == 'all'){
    if (!is.null(object$brk)){
      cbind(object$brk, object$lam)
    }else{
      message('This is a exponential model with no chage-points.')
      object$lam
    }
  }else if (parm == 'lam'){
    object$lam
  }else if (parm == 'brk'){
    object$brk
  }else {
    stop('The \'parm\' argument must be one of \'all\', \'lam\', \'brk\'.')
  }
}

summary.pwexpm <- function(object, ...){
  structure(object, class="summary.pwexpm")
}

print.summary.pwexpm <- function(x, ...){
  object <- x
  if (!is.null(object$brk)){
    s <- cbind(object$brk, object$lam)
    cat(sprintf('This is a piecewise-exponential model with %d change-point(s): \n\n', length(object$brk)))
  }else{
    s <- object$lam
    cat('This is a exponential model: \n\n')
  }
  print(cbind(s,data.frame(AIC=object$AIC, BIC=object$BIC, logLik=object$logLik)))
  cat(sprintf('The \'%s\' optimizer was used for estimation.', object$para$optimizer))
}

print.pwexpm <- function(x, ...){
  object <- x
  if (!is.null(object$brk)){
    s <- cbind(object$brk, object$lam)
    cat(sprintf('This is a piecewise-exponential model with %d change-point(s): \n\n', length(object$brk)))
  }else{
    s <- object$lam
    cat('This is a exponential model: \n\n')
  }
  print(cbind(s,data.frame(AIC=object$AIC, BIC=object$BIC, logLik=object$logLik)))
}

confint.boot.pwexpm <- function(object, parm, level = 0.90, ...){
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  if (missing(parm) || parm=='all'){
    if (!is.null(object$brk)){
      rbind(t(sapply(object$brk, function(x, level)quantile(x, probs=a))),
            t(sapply(object$lam, function(x, level)quantile(x, probs=a))))
    }else{
      message('This is a exponential model with no chage-points.')
      t(sapply(object$lam, function(x, level)quantile(x, probs=a)))
    }
  }else if (parm == 'lam'){
    t(sapply(object$lam, function(x, level)quantile(x, probs=a)))
  }else if (parm == 'brk'){
    if (!is.null(object$brk)){
      t(sapply(object$brk, function(x, level)quantile(x, probs=a)))
    }else{
      message('This is a exponential model with no chage-points.')
      NULL
    }
  }else {
    stop('The \'parm\' argument must be one of \'all\', \'lam\', \'brk\'.')
  }
}

summary.boot.pwexpm <- function(object, ...){
  structure(object, class="summary.boot.pwexpm")
}

print.summary.boot.pwexpm <- function(x, ...){
  object <- x
  if (!is.null(object$brk)){
    s <- cbind(object$brk[1,,drop=F], object$lam[1,,drop=F])
    cat(sprintf('This is a bootstrapping piecewise-exponential model with %d change-point(s): \n\n', length(object$brk)))
  }else{
    s <- object$lam[1,,drop=F]
    cat('This is a bootstrapping exponential model: \n\n')
  }
  print(cbind(s,data.frame(AIC=object$AIC, BIC=object$BIC, logLik=object$logLik)))
  cat(sprintf('The requested number of bootstrapping resampling (nsim) = %d, and %d iterations are successful. \n', object$para$nsim, NROW(object$lam)))
  cat(sprintf('The \'%s\' optimizer was used for estimation.', object$para$optimizer)) # here optimizer is from the input, may not be the real optimizer used in all iterations. Update this in future.
}


print.boot.pwexpm <- function(x, ...){
  object <- x
  if (!is.null(object$brk)){
    s <- cbind(object$brk[1,,drop=F], object$lam[1,,drop=F])
    cat(sprintf('This is a bootstrapping piecewise-exponential model with %d change-point(s): \n\n', length(object$brk)))
  }else{
    s <- object$lam[1,,drop=F]
    cat('This is a bootstrapping exponential model: \n\n')
  }
  print(cbind(s,data.frame(AIC=object$AIC, BIC=object$BIC, logLik=object$logLik)))
  cat(sprintf('The number of bootstrapping resampling (nsim) = %d.', object$para$nsim))
}

print.cv.pwexpm <- function(x, ...){
  object <- x
  cat(sprintf('The median CV log likelihood is %.3f.\n', median(object)))
  cat(sprintf('The number of resampling (nsim) = %d.', length(object)))
}

summary.predict.pwexpm <- function(object, ...){
  structure(object, class="summary.predict.pwexpm")
}

print.summary.predict.pwexpm <- function(x, ...){
  object <- x
  cat("This a predicted event curve object without bootsatraping.\n")
  cat(sprintf('The requested number of iterations (n_each) = %d, and %d simulaitons are successful. \n', object$para$n_each, length(object$event_fun)))
}

print.predict.pwexpm <- function(x, ...){
  object <- x
  cat("This a predicted event curve object without bootsatraping.\n")
  cat("Please use 'plot_event' function to plot event cruve or calculate expected events/timeline.")
}

summary.predict.boot.pwexpm <- function(object, ...){
  structure(object, class="summary.predict.boot.pwexpm")
}

print.summary.predict.boot.pwexpm <- function(x, ...){
  object <- x
  cat(sprintf('This a predicted event curve object with bootsatraping (nsim = %d).\n', object$nsim))
  cat(sprintf('The requested number of iterations for each bootstrapping sample (n_each) = %d. \n', object$para$n_each))
  cat(sprintf('The total requested number of simulations (nsim*n_each) = %d, and %d simulations are successful. \n', object$para$n_each*object$nsim, length(object$event_fun)))
}

print.predict.boot.pwexpm <- function(x, ...){
  object <- x
  cat(sprintf('This a predicted event curve object with bootsatraping (nsim = %d).\n', object$nsim))
  cat("Please use 'plot_event' function to plot event cruve or calculate expected events/timeline.")
}

AIC.pwexpm <- function(object, ..., k = 2){
  if (k!=2){
    stop('Only support the classical AIC (k=2).')
  }
  return(object$AIC)
}

BIC.pwexpm <- function(object, ...){
  return(object$BIC)
}

logLik.pwexpm <- function(object, ...){
  return(object$logLik)
}

print.sim_followup <- function(x, ...) {
  for (name in names(x)){
    if (name == 'T_all'){
      cat('Simulation result for the whole population:\n')
    }else if(name == 'T_by_group'){
      cat('Simulation result by arms:\n')
    }else if(name == 'T_by_strata'){
      cat('Simulation result by strata:\n')
    }else if(name == 'T_by_group_strata'){
      cat('Simulation result by arms and strata:\n')
    }
    print(x[[name]])
    cat('\n')
  }
}
