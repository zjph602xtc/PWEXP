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

print.predict.pwexpm <- function(x, ...){
  object <- x
  cat("This a predicted event curve object without bootsatraping.\n")
  cat("Please use 'plot_event' function to plot event cruve or calculate expected events/timeline.")
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
