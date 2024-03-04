rpwexp <- function(n, rate=1, breakpoint=NULL){
  x <- runif(n)
  return(PWEXP::qpwexp(x, rate, breakpoint, lower.tail=TRUE, log.p=FALSE))
}

rpwexp_conditional <- function(n, qT, rate, breakpoint=NULL){
  # length of qT should be 1 or == length(n)
  x <- runif(n)
  return(qpwexp_conditional(x, qT, rate, breakpoint, lower.tail=TRUE, log.p=FALSE))
}

ppwexp <- function(q, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # S(t)=exp(-((lam1-lam2)*d1+(lam2-lam3)*d2+... lami*t))
  # breakpoint must be sorted!
  # if safety_check is FALSE, then skip all checks and one_piece must be correctly specified.
  if (safety_check){
    if (!is.null(breakpoint) & length(breakpoint)==0){
      stop('breakpoint must be NULL or a vector of numerics.')
    }
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
  }
  if (one_piece){
    return(pexp(q, rate, lower.tail, log.p))
  }else {
    shift_t <- c(0, cumsum(-diff(rate)*breakpoint))
    interval <- findInterval(q, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
    logs <- -(rate[interval]*q+shift_t[interval])
    if (!lower.tail & log.p){
      return(logs)
    }else if (!lower.tail & !log.p){
      return(exp(logs))
    }else if (lower.tail & !log.p){
      return(1-exp(logs))
    }else {
      return(log(1-exp(logs)))
    }
  }
}

ppwexp_conditional <- function(q, qT, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # F(t|X>qT) = 1 - S(t)/S(qT)
  # S(t|X>qT) = S(t)/S(qT)
  # breakpoint must be sorted!
  # if safety_check is FALSE, then skip all checks and one_piece must be correctly specified.
  if (safety_check){
    if (!is.null(breakpoint) & length(breakpoint)==0){
      stop('breakpoint must be NULL or a vector of numerics.')
    }
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
    if (length(qT)!=1 & length(q)!=length(qT)){
      stop('length of qT should be 1 or == length(q)')
    }
    if (any(q < qT)){
      stop('q must be equal or larger than qT')
    }
  }
  logs <- PWEXP::ppwexp(q, rate, breakpoint, lower.tail = F, log.p = T, one_piece, safety_check = FALSE)-
    PWEXP::ppwexp(qT, rate, breakpoint, lower.tail = F, log.p = T, one_piece, safety_check = FALSE)
  if (!lower.tail & log.p){
    return(logs)
  }else if (!lower.tail & !log.p){
    return(exp(logs))
  }else if (lower.tail & !log.p){
    return(1-exp(logs))
  }else {
    return(log(1-exp(logs)))
  }
}

dpwexp <- function(x, rate=1, breakpoint=NULL, log=FALSE, one_piece, safety_check=TRUE){
  # f(t)=lambda_i * exp(-((lam1-lam2)*d1+(lam2-lam3)*d2+... lami*t))
  # breakpoint must be sorted!
  if (safety_check){
    if (!is.null(breakpoint) & length(breakpoint)==0){
      stop('breakpoint must be NULL or a vector of numerics.')
    }
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
  }
  if (one_piece){
    return(dexp(x, rate, log))
  }else {
    shift_t <- c(0, cumsum(-diff(rate)*breakpoint))
    log_rate <- log(rate)
    interval <- findInterval(x, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
    logd <- -rate[interval]*x-(shift_t-log_rate)[interval]
    # logd <- -(rate[interval]*x+shift_t[interval])+log(rate[interval])
    if (log){
      return(logd)
    }else {
      return(exp(logd))
    }
  }
}

qpwexp <- function(p, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # F(y)^(-1)=(-log(1-y)-(lam1-lam2)*d1-(lam2-lam3)*d2-...)/lami
  # breakpoint must be sorted!
  if (safety_check){
    if (!is.null(breakpoint) & length(breakpoint)==0){
      stop('breakpoint must be NULL or a vector of numerics.')
    }
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
  }
  if (one_piece){
    return(qexp(p, rate, lower.tail, log.p))
  }else {
    if (log.p){
      p <- exp(p)
    }
    shift_t <- c(0, cumsum(-diff(rate)*breakpoint))
    Fcut <- PWEXP::ppwexp(breakpoint, rate, breakpoint, lower.tail, log.p=FALSE, one_piece, safety_check = FALSE)
    if (lower.tail){
      interval <- findInterval(p, vec=c(-Inf, Fcut, Inf), rightmost.closed = FALSE)
      t <- (-log(1-p)-shift_t[interval])/rate[interval]
    }else{
      interval <- 2+length(Fcut)-findInterval(p, vec=rev(c(Inf, Fcut, -Inf)), left.open = TRUE, rightmost.closed = TRUE)
      t <- (-log(p)-shift_t[interval])/rate[interval]
    }
    return(t)
  }
}

qpwexp_conditional <- function(p, qT, rate=1, breakpoint=NULL, lower.tail=TRUE, log.p=FALSE, one_piece, safety_check=TRUE){
  # t=F(y)^(-1) given X>=qT
  # assume qT in Mth interval, t (or say p) in kth interval
  # t=(-log(1-y)-(lam_{M}-lam_{M+1})*d_M-...-(lam_{k-1}-lam_{k})*d_{k-1}+lam_M*T)/lam_{K}
  # breakpoint must be sorted!

  # length of qT should be 1 or == length(p)
  if (safety_check){
    if (!is.null(breakpoint) & length(breakpoint)==0){
      stop('breakpoint must be NULL or a vector of numerics.')
    }
    one_piece <- is.null(breakpoint)
    if ((length(rate)-1) != length(breakpoint)){
      stop('wrong number of breakpoints')
    }
    if (any(rate < 0)){
      stop('negative rate')
    }
    if (any(diff(breakpoint)<=0)) {
      stop('breakpoint must be sorted')
    }
    if (length(qT)!=1 & length(p)!=length(qT)){
      stop('length of qT should be 1 or == length(p)')
    }
  }
  if (log.p){
    p <- exp(p)
  }
  if (one_piece){
    if (lower.tail){
      t <- (-log(1-p)+rate*qT)/rate
    }else {
      t <- (-log(p)+rate*qT)/rate
    }
    return(t)
  }
  if (length(qT)!=1){
    res <- mapply(function(p,qT)qpwexp_conditional(p, qT, rate, breakpoint, lower.tail, log.p, one_piece, safety_check = FALSE), p=p,qT=qT)
    return(res)
  }else {
    shift_t <- -diff(rate)*breakpoint
    Fcut <- ppwexp_conditional(breakpoint, qT, rate, breakpoint, lower.tail, log.p=FALSE, one_piece, safety_check = FALSE)
    if (lower.tail){
      interval_k <- findInterval(p, vec=c(-Inf, Fcut, Inf), rightmost.closed = FALSE)
      interval_m <- findInterval(qT, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
      # shift_int <- mapply(function(i,j)shift_t[i:j], i=interval_m, j=interval_k-1, SIMPLIFY = FALSE)
      shift_int <- lapply(interval_k-1, function(k)shift_t[interval_m:k])
      shift_int[interval_m > (interval_k-1)] <- 0
      shift_int <- sapply(shift_int, sum)
      t <- (-log(1-p)-shift_int+rate[interval_m]*qT)/rate[interval_k]
    }else{
      interval_k <- 2+length(Fcut)-findInterval(p, vec=rev(c(Inf, Fcut, -Inf)), left.open = TRUE, rightmost.closed = TRUE)
      interval_m <- findInterval(qT, vec=c(-Inf, breakpoint, Inf), rightmost.closed = FALSE)
      # shift_int <- mapply(function(i,j)shift_t[i:j], i=interval_m, j=interval_k-1, SIMPLIFY = FALSE)
      shift_int <- lapply(interval_k-1, function(k)shift_t[interval_m:k])
      shift_int[interval_m > (interval_k-1)] <- 0
      shift_int <- sapply(shift_int, sum)
      t <- (-log(p)-shift_int+rate[interval_m]*qT)/rate[interval_k]
    }
    return(t)
  }
}
