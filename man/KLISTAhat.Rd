\name{KLISTAhat}
\alias{KLISTAhat}

\title{Estimation of the Space-Time Inhomogeneous K LISTA functions}

\description{Compute an estimate of the space-time K LISTA functions.}

\usage{
KLISTAhat(xyt, s.region, t.region, dist, times, lambda, correction = "isotropic")
}

\arguments{
  \item{xyt}{Coordinates and times \eqn{(x,y,t)}{(x,y,t)} of the point pattern.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  
  \item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  
  \item{dist}{Vector of distances \eqn{u}{u} at which \eqn{K^{(i)}(u,v)}{K^{(i)}(u,v)} is computed. If missing, the maximum of \code{dist} is given by \eqn{\min(S_x,S_y)/4}{min(S_x, S_y)/4}, where \eqn{S_x}{S_x} and \eqn{S_y}{S_y} represent the maximum width and height of the bounding box of \code{s.region}.}
  
  \item{times}{Vector of times \eqn{v}{v} at which \eqn{K^{(i)}(u,v)}{K^{(i)}(u,v)} is computed.  If missing, the maximum of \code{times} is given by \eqn{(T_{\max} - T_{\min})/4}{(T_max - T_min)/4}, where \eqn{T_{\min}}{T_min} and \eqn{T_{\max}}{T_max} are the minimum and maximum of the time interval \eqn{T}{T}.}
  
  \item{lambda}{Vector of values of the space-time intensity function evaluated at the points \eqn{(x,y,t)}{(x,y,t)} in \eqn{S\times T}{S x T}. If \code{lambda} is missing, the estimate of the space-time pair correlation function is computed as for the homogeneous case, i.e. considering \eqn{(n-1)/|S \times T|}{(n-1)/|S x T|} as an estimate of the space-time intensity.}
  

  \item{correction}{A character vector specifying the edge correction(s) to be applied among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"} and \code{"none"} (see \code{\link{PCFhat}}). The default is \code{"isotropic"}.
}
}

\details{An individual product density LISTA functions \eqn{K^{(i)}(.,.)}{K^{(i)}(.,.)} should reveal the extent of the contribution of the event \eqn{(u_i,t_i)}{(u_i,t_i)} to the global estimator of the K-function \eqn{K(.,.)}{K(.,.)}, and may provide a further description of structure in the data (e.g., determining events with similar local structure through dissimilarity measures of the individual LISTA functions), for more details see Siino et al. (2019).}

\value{
A list containing:

\item{list.KLISTA}{ A list containing the values of the estimation of \eqn{K^{(i)}(r,t)}{K^{(i)}(r,t)} for each one of \eqn{n} points of the point pattern by matrixs.}
  
\item{klistatheo}{\code{ndist} x \code{ntimes} matrix containing theoretical values for a Poisson process.}
  
\item{dist, times}{Parameters passed in argument.}
   
\item{correction}{The name(s) of the edge correction method(s) passed in argument.}

}

\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.

Siino, M., Adelfio, G., Mateu, J. and Rodriguez-Cortes, F. J. (2019). Some properties of weighted local second-order statistcs for spatio-temporal point process. Submitted.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
}
\author{Francisco J. Rodriguez-Cortes <frrodriguezc@unal.edu.co>}
