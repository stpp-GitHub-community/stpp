\name{plot.stpp}
\alias{plot.stpp}
\alias{plot}

\title{Plot for spatio-temporal point objects}

\description{
This function plot either xy-locations and cumulative distribution of the times
or xy-locations with time treated as a quantitative mark attached to each location.}

\usage{
\method{plot}{stpp}(x, s.region=NULL, t.region=NULL, mark=FALSE, mark.cexmin=0.4,
     mark.cexmax=1.2, mark.col=1, ...)
} 

\arguments{
  \item{x}{any object of class \code{stpp} in spatio-temporal point format.}
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations. If \code{s.region} is missing, the default limits are
  considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing,  the default limits are
  considered.}
  \item{mark}{Logical. If \code{FALSE} (default), xy-locations and cumulative 
  distribution of the times are plotted. If \code{TRUE}, the time is treated as 
  a quantitative mark attached to each location, and the locations are plotted 
  with the size and/or colour of the plotting symbol determined by the value of 
  the mark.}
  \item{mark.cexmin, mark.cexmax}{range of the size of the plotting symbol when
  \code{mark=TRUE}.}
  \item{mark.col}{colour of the plotting symbol when \code{mark=TRUE}. If \code{mark.col=0}, 
  all locations have the same colour specified by the usual \code{col} argument.
  Otherwise, can be 1 or "black" (default), 2 or "red", 3 or "green", 4 or "blue", in which cases
  symbols colour is faded, and the darker corresponds to the most recent time.}
  \item{...}{further arguments to be passed to the function plot.
  Typical argument is \code{pch}.}
 }

\value{
None}

\seealso{
 \code{\link{as.3dpoints}} for creating data in spatio-temporal point format.
 }

\author{Edith Gabriel <edith.gabriel@univ-avignon.fr>.}

\references{
Gabriel E., Rowlingson B., Diggle P. (2013)
stpp: an R package for plotting, simulating and analyzing 
Spatio-Temporal Point Patterns. Journal of Statistical Software, 
53(2), 1--29.
}
