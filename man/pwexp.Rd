\name{piecewise exponential}
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
  \item{breakpoint}{a vector of breakpoints. The length is \code{length(rate)-1}. Can be NULL if rate is a single value.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{one_piece}{(only required when \code{safety_check=FALSE}) whether the distribution only has one piece (i.e.,
  rate is a single value and breakpoint=NULL).}
  \item{safety_check}{logical; whether check the input arguments; if FALSE, function has better computing performance by skipping all safety checks.}
}

\value{
\code{dpwexp} gives the density, \code{ppwexp} gives the distribution function, \code{qpwexp} gives the quantile function, and \code{rpwexp} generates random deviates.

The length of the result is determined by \code{x}, \code{q}, \code{p} or \code{n} for \code{dpwexp}, \code{ppwexp}, \code{qpwexp} or \code{rpwexp}. You can only specify a single piecewise exponential distribution every time you call these functions. This is different from the exponential distribution functions in package \pkg{stats}.

Arguments \code{rate} and \code{breakpoint} must match. The length of rate is the length of breakpoint + 1.
}
\details{
The piecewise distribution function with piecewise rate \eqn{\lambda_1, \dots, \lambda_r} is
\deqn{f(t)=\lambda_{r+1} exp{[\sum_{i=1}^r(\lambda_{i+1}-\lambda_{i})d_i-\lambda_{r+1}t}}
for \eqn{d_{r} \le t < d_{r+1}}.

See webpage \url{https://zjph602xtc.github.io/PWEXP/} for more details for its hazard function, cumulative hazard function, survival function, cumulative density function, quantile function.
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
# use rpwexp function to generate piecewise exp samples with rate 2, 1, 3
r_sample <- rpwexp(50000, rate=c(2, 1, 3), breakpoint=c(0.3, 0.8))
hist(r_sample, freq=FALSE, breaks=200, main="Density of Piecewsie Exp Dist",
    xlab='t', xlim=c(0, 1.2))

# piecewise exp density with rate 2, 1, 3
t <- seq(0, 1.5, 0.01)
f2 <- dpwexp(t, rate=c(2, 1, 3), breakpoint=c(0.3, 0.8))
points(t, f2, col='red', pch=16)

# exp distribution can be a special case of piecewise exp distribution
f1 <- dpwexp(t, rate=2)
lines(t, f1, lwd=2)
legend('topright', c('exp dist with rate 2','piecewise exp dist with rate 2, 1,
    3','histogram of piecewise exp dist with rate 2, 1, 3'),
    col=c('black','red'), fill=c(NA, NA, 'grey'), border=c('white', 'white',
    'black'), lty=c(1, NA, NA), pch=c(NA, 16, NA), lwd=2)

# CDF of piecewise exp with rate 2, 1, 3
F2 <- ppwexp(t, rate=c(2, 1, 3), breakpoint=c(0.3, 0.8), lower.tail=TRUE)
plot(t, F2, type='l', col='red', lwd=2, main="CDF and Quantile Function of
    Piecewsie Exp Dist", xlim=c(0, 1.5), ylim=c(0, 1.5))

# CDF of exp dist is compatible with our package
F1 <- ppwexp(t, rate=2, lower.tail=TRUE)
lines(t, F1, lwd=2)

# plot quantile functions of both distributions
lines(F1, qpwexp(F1, rate=2, lower.tail=TRUE), lty=2, lwd=2)
lines(F2, qpwexp(F2, rate=c(2, 1, 3), breakpoint=c(0.3,0.8), lower.tail=TRUE),
    col='red', lty=2, lwd=2)

abline(0, 1, col='grey')
legend('topleft', c('CDF of piecewise exp with rate 2, 1, 3', 'quantile
    function of piecewise exp with rate 2, 1, 3', 'CDF of exp with rate 2',
    'quantile function of exp with rate 2'), col=c('red', 'red', 'black',
    'black'), lty=c(1, 2, 1, 2), lwd=2)
}










