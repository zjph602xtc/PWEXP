\name{confint.boot.pwexpm}
\alias{confint.boot.pwexpm}
\title{Confidence Intervals for Model Parameters in boot.pwexpm}
\description{
  Computes confidence intervals for one or more parameters in a bootstrapped PWE model.
}

\usage{
\method{confint}{boot.pwexpm}(object, parm, level = 0.90, ...)
}
\arguments{
 \item{object}{a \code{boot.pwexpm} object.}
 \item{parm}{a specification of which parameters are to be given confidence intervals (either \code{'lam'} for hazard rates, or \code{'brk'} for change-points). If missing, all parameters are considered.}
 \item{level}{the confidence level required.}
 \item{...}{some methods for this generic require additional arguments. None are used in this method.}
}

\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\value{
 A matrix with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in \%.
}

