\name{plot.PICBayes}
\alias{plot.PICBayes}
\title{Plot method for a PICBayes model}
\description{Plot estimated baseline survival function at \code{grids}.}
\usage{
\method{plot}{PICBayes}(x, y, ...)}
\arguments{
  \item{x}{A sequence of points (\code{grids}) where baseline survival probabilities 
	are estimated.}
  \item{y}{Estiamted baseline survival at \code{grids}.}
  \item{...}{Other arguments if any.}
}
\value{A plot of baseline survival function.}
