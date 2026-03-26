\name{coef.pwexpm}
\alias{coef.pwexpm}
\title{Extract pwexpm Coefficients}
\description{
  The estimated coefficients for the \code{pwexpm} are extracted.
}
\usage{
\method{coef}{pwexpm}(object, parm, ...)
}
\arguments{
 \item{object}{a \code{pwexpm} object.}
 \item{parm}{a specification of which parameters are to be extracted (either \code{'lam'} for hazard rates, or \code{'brk'} for change-points). If missing, all parameters are considered.}
 \item{...}{some methods for this generic require additional arguments. None are used in this method.}
}

\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\value{
 A data frame containing the requested coefficients.
}

