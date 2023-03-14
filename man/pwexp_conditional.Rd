\name{conditional piecewise exponential}
\alias{ppwexp_conditional}
\alias{qpwexp_conditional}
\alias{rpwexp_conditional}
\title{The Conditional Piecewise Exponential Distribution}
\description{
Distribution function, quantile function and random generation for the piecewise exponential distribution \eqn{t} with piecewise rate \code{rate} given \eqn{t>qT}.
}
\usage{
ppwexp_conditional(q, qT, rate=1, breakpoint=NULL, lower.tail=TRUE,
                   log.p=FALSE, one_piece, safety_check=TRUE)
qpwexp_conditional(p, qT, rate=1, breakpoint=NULL, lower.tail=TRUE,
                   log.p=FALSE, one_piece, safety_check=TRUE)
rpwexp_conditional(n, qT, rate, breakpoint=NULL)
}
\arguments{
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{qT}{the distribution is conditional on \eqn{t>}\code{qT}. \code{qT} can be a scalar or a vector with the same length of \code{q} or \code{p}.}
  \item{n}{number of observations. Must be a positive integer with length 1. }
  \item{rate}{a vector of rates in each piece.}
  \item{breakpiont}{a vector of breakpoints. The length is \code{length(rate)-1}. Can be NULL if rate is a single value.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{one_piece}{(only required when \code{safety_check=FALSE}) whether the distribution only has one piece (i.e.,
  rate is a single value and breakpoint=NULL).}
  \item{safety_check}{logical; whether check the input arguments; if FALSE, function has better computing performance by skipping all safety checks.}
}

\value{
\code{ppwexp_conditional} gives the conditional distribution function, \code{qpwexp_conditional} gives the conditional quantile function, and \code{rpwexp_conditional} generates conditional random variables.

The length of the result is determined by \code{q}, \code{p} or \code{n} for \code{ppwexp_conditional}, \code{qpwexp_conditional} or \code{rpwexp_conditional}. You can only specify a single piecewise exponential distribution every time you call these functions. This is different from the exponential distribution functions in package \pkg{stats}.

When the length of \code{qT} is 1, then all results are conditional on the same \eqn{t>}\code{qT}.
In \code{rpwexp_conditional}, \code{qT} must be a scalar. When the length of \code{qT} equals to the length of \code{q} or \code{p}, then each value in the result is conditional on \eqn{t>}\code{qT} for each value in \code{qT}.

Arguments \code{rate} and \code{breakpoint} must match. The length of rate is the length of breakpoint + 1.
}
\details{
See webpage \url{https://zjph602xtc.github.io/PWEXP/} for more details for its survival function, cumulative density function, quantile function.
}


\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{
  \code{\link{dpwexp}},
  \code{\link{ppwexp}},
  \code{\link{qpwexp}},
  \code{\link{rpwexp}}
}
\examples{
t <- seq(0,1.2,0.01)

# distribution function
F2 <- ppwexp(t, rate = c(2, 1, 3), breakpoint = c(0.3,0.8), lower.tail = F)
plot(t, F2, type='l', lwd=2, xlim=c(0,1.2), ylim=c(0,1.2))
t <- seq(0.1,1.2,0.01)
F2_con <- ppwexp_conditional(t, qT = 0.1, rate = c(2, 1, 3), breakpoint = c(0.3,0.8),
                             lower.tail = F)
lines(t, F2_con, col='red', lwd=2)
abline(v=0.1,col='grey')

# quantile function
abline(0,1,col='grey')
lines(F2, qpwexp(F2, rate = c(2, 1, 3), breakpoint = c(0.3,0.8), lower.tail = F), lty=2, lwd=2)
lines(F2_con, qpwexp_conditional(F2_con, qT = 0.1, rate = c(2, 1, 3), breakpoint = c(0.3,0.8),
                                 lower.tail = F), lty=2, lwd=2, col='red')

# data generation
plot(t, 1-F2_con, type='l', lwd=2, xlim=c(0,1.2), ylim=c(0,1))
lines(ecdf(rpwexp_conditional(3000, qT = 0.1, rate = c(2, 1, 3),
                              breakpoint = c(0.3,0.8))), col='red', lwd=2)

# when qT is a vector
t <- seq(0.1,1.2,0.1)
ppwexp_conditional(t, qT = seq(0.1,0.5,length=12), rate = c(2, 1, 3),
                   breakpoint = c(0.3,0.8), lower.tail = F)
}










