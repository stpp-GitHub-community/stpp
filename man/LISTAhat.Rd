\name{LISTAhat}
\alias{LISTAhat}

\title{Estimation of the Space-Time Inhomogeneous Pair Correlation LISTA functions}

\description{Compute an estimate of the space-time pair correlarion LISTA functions.}

\usage{
LISTAhat(xyt, s.region, t.region, dist, times, lambda,
ks = "box", hs, kt = "box", ht, correction = "isotropic")
}

\arguments{
  \item{xyt}{Coordinates and times \eqn{(x,y,t)}{(x,y,t)} of the point pattern.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  
  \item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  
  \item{dist}{Vector of distances \eqn{u}{u} at which \eqn{g^{(i)}(u,v)}{g^{(i)}(u,v)} is computed. If missing, the maximum of \code{dist} is given by \eqn{\min(S_x,S_y)/4}{min(S_x, S_y)/4}, where \eqn{S_x}{S_x} and \eqn{S_y}{S_y} represent the maximum width and height of the bounding box of \code{s.region}.}
  
  \item{times}{Vector of times \eqn{v}{v} at which \eqn{g^{(i)}(u,v)}{g^{(i)}(u,v)} is computed.  If missing, the maximum of \code{times} is given by \eqn{(T_{\max} - T_{\min})/4}{(T_max - T_min)/4}, where \eqn{T_{\min}}{T_min} and \eqn{T_{\max}}{T_max} are the minimum and maximum of the time interval \eqn{T}{T}.}
  
  \item{lambda}{Vector of values of the space-time intensity function evaluated at the points \eqn{(x,y,t)}{(x,y,t)} in \eqn{S\times T}{S x T}. If \code{lambda} is missing, the estimate of the space-time pair correlation function is computed as for the homogeneous case, i.e. considering \eqn{(n-1)/|S \times T|}{(n-1)/|S x T|} as an estimate of the space-time intensity.}
  
  \item{ks}{Kernel function for the spatial distances. Default is the \code{"box"} kernel. Can also be \code{"epanech"} for the Epanechnikov kernel or \code{"gaussian"} or \code{"biweight"}.}
  
  \item{hs}{Bandwidth of the kernel function \code{ks}.}
  
  \item{kt}{Kernel function for the temporal distances. Default is the \code{"box"} kernel. Can also be \code{"epanech"} for the Epanechnikov kernel or \code{"gaussian"} or \code{"biweight"}.}
  
  \item{ht}{Bandwidth of the kernel function \code{kt}.}
  
  \item{correction}{A character vector specifying the edge correction(s) to be applied among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"} and \code{"none"} (see \code{\link{PCFhat}}). The default is \code{"isotropic"}.
}
}

\details{An individual product density LISTA functions \eqn{g^{(i)}(.,.)}{g^{(i)}(.,.)} should reveal the extent of the contribution of the event \eqn{(u_i,t_i)}{(u_i,t_i)} to the global estimator of the pair correlation function \eqn{g(.,.)}{g(.,.)}, and may provide a further description of structure in the data (e.g., determining events with similar local structure through dissimilarity measures of the individual LISTA functions), for more details see Siino et al. (2017).}

\value{
A list containing:

\item{list.LISTA}{ A list containing the values of the estimation of \eqn{g^{(i)}(r,t)}{g^{(i)}(r,t)} for each one of \eqn{n} points of the point pattern by matrixs.}
  
\item{listatheo}{\code{ndist} x \code{ntimes} matrix containing theoretical values for a Poisson process.}
  
\item{dist, times}{Parameters passed in argument.}

\item{kernel}{Vector of names and bandwidths of the spatial and temporal kernels.}
   
\item{correction}{The name(s) of the edge correction method(s) passed in argument.}

}

\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.

Siino, M., Rodriguez-Cortes, F. J., Mateu, J. and Adelfio, G. (2017). Testing for local structure in spatio-temporal point pattern data. Environmetrics. DOI: 10.1002/env.2463.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
}
\author{Francisco J. Rodriguez-Cortes <frrodriguezc@unal.edu.co>}
