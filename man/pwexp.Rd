\name{Piecewise exponential}
\alias{dpwexp}
\alias{ppwexp}
\alias{qpwexp}
\alias{rpwexp}
\title{The Piecewise Exponential Distribution}
\description{
Density, distribution function, quantile function and random generation for the piecewise exponential distribution with piecewise rate \code{rate}.
}
\usage{
dpwexp(x, rate = 1, breakpoint = NULL, log = FALSE, one_piece, safety_check = TRUE)
ppwexp(q, rate = 1, breakpoint = NULL, lower.tail = TRUE, log.p = FALSE,
       one_piece, safety_check = TRUE)
qpwexp(p, rate = 1, breakpoint = NULL, lower.tail = TRUE, log.p = FALSE,
       one_piece, safety_check = TRUE)
rpwexp(n, rate = 1, breakpoint = NULL)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
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
\code{dpwexp} gives the density, \code{ppwexp} gives the distribution function, \code{qpwexp} gives the quantile function, and \code{rpwexp} generates random deviates.

The length of the result is determined by \code{x}, \code{q}, \code{p} or \code{n} for \code{dpwexp}, \code{ppwexp}, \code{qpwexp} or \code{rpwexp}. This is different from the exponential distribution functions in package \pkg{stats}.

Arguments \code{rate} and \code{breakpoint} must match, and only determine a single piecewise distribution.
}
\details{
dd
}

\note{
When \code{breakpoint=NULL}, the function calls exponential function in \pkg{stats}.
}

\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\seealso{
  \code{\link{ppwexp_conditional}},
  \code{\link{qpwexp_conditional}},
  \code{\link{rpwexp_conditional}}
}
\examples{
t <- seq(0,1.2,0.01)

# density function
f1 <- dpwexp(t, rate = 2)
plot(t, f1, type='l', lwd=2)
f2 <- dpwexp(t, rate = c(2, 1, 3), breakpoint = c(0.3,0.8))
points(t, f2, col='red', pch=16)

# generation function
hist(rpwexp(50000, rate = c(2, 1, 3), breakpoint = c(0.3,0.8)), freq = F, add = T, breaks = 500)

# distribution function
F1 <- ppwexp(t, rate = 2, lower.tail = F)
plot(t, F1, type='l', lwd=2, xlim=c(0,1.2), ylim=c(0,1.2))
F2 <- ppwexp(t, rate = c(2, 1, 3), breakpoint = c(0.3,0.8), lower.tail = F)
lines(t, F2, col='red', lwd=2)

# quantile function
abline(0,1,col='grey')
lines(F1, qpwexp(F1, rate = 2, lower.tail = F), lty=2, lwd=2)
lines(F2, qpwexp(F2, rate = c(2, 1, 3), breakpoint = c(0.3,0.8), lower.tail = F), col='red', lty=2, lwd=2)
}










