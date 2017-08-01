\name{sim.stpp}
\alias{sim.stpp}
\title{Generate spatio-temporal point patterns}

\description{
Generate one (or several) realisation(s) of a spatio-temporal point process in a region S x T.
}

\usage{
sim.stpp(class="poisson", s.region, t.region, npoints=NULL, 
nsim=1, ...)
}

\arguments{
  \item{class}{Must be chosen among "poisson", "cluster", "cox", "infectious" and "inhibition".} 
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations. 
   If \code{s.region} is missing, the unit square is considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. 
   If \code{t.region} is missing, the interval [0,1] is considered.}
  \item{npoints}{number of points to simulate.} 
  \item{nsim}{number of simulations to generate. Default is 1.}
  \item{...}{additional parameters related to the \code{class} parameter.
  See \code{\link{rpp}} for the Poisson process; \code{\link{rpcp}} for the Poisson 
  cluster process; \code{\link{rlgcp}} for the Log-Gaussian Cox process; 
  \code{\link{rinter}} for the interaction (inhibition or contagious) 
  process and \code{\link{rinfec}} for the infectious process.}
}

\value{
A list containing:
\item{xyt}{matrix (or list of matrices if \code{nsim}>1)
containing the points (x,y,t) of the simulated point pattern.
\code{xyt} (or any element of the list if \code{nsim}>1) is an object 
of the class \code{stpp}.}
\item{s.region, t.region}{parameters passed in argument.}
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr>
}

\seealso{
 \code{\link{rpp}}, \code{\link{rinfec}}, \code{\link{rinter}}, 
 \code{\link{rpcp}} and  \code{\link{rlgcp}} for the simulation of 
 Poisson, infectious, interaction, Poisson cluster and log-gaussian 
 Cox processes respectively; and \code{\link{plot.stpp}}, 
 \code{\link{animation}} and \code{\link{stan}} for plotting space-time 
 point patterns. }

