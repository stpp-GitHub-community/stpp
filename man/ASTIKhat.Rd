\name{ASTIKhat}
\alias{ASTIKhat}
\title{Anisotropic space-time inhomogeneous \eqn{K}-function}
\usage{
ASTIKhat(xyt, s.region, t.region, lambda, dist, times, ang,
  correction = "border")
}
\arguments{
\item{xyt}{Coordinates and times \eqn{(x,y,t)} of the point pattern.}

\item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If s.region is missing, the Ripley-Rasson estimate of the spatial domain is considered.}

\item{t.region}{Vector containing the minimum and maximum values of the time interval. If t.region is missing, the range of \eqn{xyt[,3]} is considered.}

\item{lambda}{Vector of values of the spatio-temporal intensity function evaluated at the points \eqn{(x,y,t)} in \eqn{W x T}. If lambda is missing, the estimate of the anisotropic space-time \eqn{K}-function is computed as for the homogeneous case, i.e. considering \eqn{n/|W x T|} as an estimate of the space-time intensity.}

\item{dist}{Vector of distances \eqn{u} at which \eqn{\hat{K}_{\phi}(r,t)} is computed.}

\item{times}{Vector of times \eqn{v} at which \eqn{\hat{K}_{\phi}(r,t)} is computed.}

\item{ang}{Angle in radians at which \eqn{\hat{K}_{\phi}(r,t)} is computed.}

\item{correction}{A character vector specifying the edge correction(s) to be applied among "border", "modified.border", "translate" and "none". The default is "border".}
}
\value{
A list containing:
\itemize{
\item \code{AKhat}: Matrix contain \code{ndist} by \code{ntimes} elements of \eqn{\hat{K}_{\phi}(r,t)}.
\item \code{dist}: If \code{dist} is missing, a vector of distances \code{u} at which \eqn{\hat{K}_{\phi}(u,v)} is computed by default.
\item \code{times}: If \code{times} is missing, a vector of distances \code{v} at which \eqn{\hat{K}_{\phi}(u,v)} is computed by default.
\item \code{Correction}: Value passed in the initial parameters.
}
}
\description{
This function compute an estimation of the anisotropic space-time inhomogeneous \eqn{K}-function.
}
\examples{
## Not run:
#################

# Realisations of the homogeneous space-time Poisson processes
stp <- rpp(100)$xyt
# Generated space-time point pattern
plot(stp)

out <- ASTIKhat(stp,ang=pi/6)

#################
}
\author{
Francisco J. Rodriguez-Cortes <cortesf@uji.es> 
}
\references{
Gabriel, E., Rowlingson, B., Diggle P J. (2013). \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.

Illian, J. B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Gonzalez, J. A., Rodriguez-Cortes, F. J., Cronie, O., Mateu, J. (2016). Spatio-temporal point process statistics: a review. Spatial Statistics. Accepted.

Ohser, J. and D. Stoyan (1981). On the second-order and orientation analysis of planar stationary point processes. Biometrical Journal 23, 523-533.
}