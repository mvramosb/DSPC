% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{parcoord2}
\alias{parcoord2}
\title{Parallel Coordinates Plot}
\usage{
parcoord2(
  x,
  col = 1,
  lty = 1,
  var.label = FALSE,
  liminf = NULL,
  limsup = NULL,
  ...
)
}
\arguments{
\item{x}{A multivariate data matrix representing reference batches or
a testing batch.}

\item{col}{A vector of colours, recycled as necessary for each observation.}

\item{lty}{A vector of line types, recycled as necessary for each observation.}

\item{var.label}{If TRUE, each variable's axis is labelled with maximum and minimum values.}

\item{liminf}{a numeric vector conformed by the inferior limit for every variable.}

\item{limsup}{a numeric vector conformed by the inferior limit for every variable.}

\item{...}{Further graphics parameters which are passed to \code{matplot}}
}
\description{
A modified version of the Parallel Coordinates plotting function
available in the MASS package. If number of rows in the table
is greater than 500, a sample is picked.
}
\details{
This function creates a parallel coordinates plot using the \code{matplot}
function. No assignable output is created.
}
\examples{
data( PlasticBags )
x = PlasticBags$Add$data[1:100,]
colorset = c(rep("black",50),rep("blue",50))
Xbar = colMeans(PlasticBags$Ref$data); Desv = apply(PlasticBags$Ref$data,2,sd)
inferiorlim = Xbar - 5*Desv
superiorlim = Xbar + 5*Desv
parcoord2( x, var.label = TRUE , col= colorset ,liminf = inferiorlim, limsup = superiorlim  )
}
\references{
Wegman, E. J. (1990) Hyperdimensional data analysis using parallel
coordinates. Journal of the American Statistical Association 85, 664 - 675.

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S.
Fourth edition. Springer.
}
\author{
B. D. Ripley. Enhancements based on ideas and code by Fabian Scheipl.
}
