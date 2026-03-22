\name{coef.pwexpm}
\alias{coef.pwexpm}
\title{Extract pwexpm Coefficients}
\description{
  The estimated coefficients for the \code{pwexpm} are extracted.
}
\usage{
\method{coef}{pwexpm}(object, parm = 'all', ...)
}
\arguments{
 \item{object}{a \code{pwexpm} object.}
 \item{parm}{a specification of which parameters are to be extracted: \code{all} for both change-points and hazard rates (lambda), \code{lam} for hazard rates (lambda), or \code{brk} for change-points.}
 \item{...}{some methods for this generic require additional arguments. None are used in this method.}
}

\author{
Tianchen Xu \email{zjph602xutianchen@gmail.com}
}
\value{
 A data frame containing the requested coefficients.
}

