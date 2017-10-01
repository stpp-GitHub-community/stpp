\name{ASTIKhat}
\alias{ASTIKhat}
\title{Anisotropic spatio-temporal inhomogeneous \eqn{K}-function}
\usage{
ASTIKhat(xyt, s.region, t.region, lambda, dist, times, ang,
  correction = "border")
}
\arguments{
\item{xyt}{Coordinates and times \eqn{(x,y,t)} of the point pattern.}

\item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If s.region is missing, the Ripley-Rasson estimate of the spatial domain is considered.}

\item{t.region}{Vector containing the minimum and maximum values of the time interval. If t.region is missing, the range of \eqn{xyt[,3]} is considered.}

\item{lambda}{Vector of values of the spatio-temporal intensity function evaluated at the points \eqn{(x,y,t)} in \eqn{W x T}. If lambda is missing, the estimate of the anisotropic spatio-temporal \eqn{K}-function is computed as for the homogeneous case, i.e. considering \eqn{n/|W x T|} as an estimate of the spatio-temporal intensity.}

\item{dist}{Vector of distances \eqn{u} at which \eqn{\hat{K}_{\phi}(r,t)} is computed.}

\item{times}{Vector of times \eqn{v} at which \eqn{\hat{K}_{\phi}(r,t)} is computed.}

\item{ang}{Vector of angles in radians at which \eqn{\hat{K}_{\phi}(r,t)} is computed.  If ang is missing, the function HASTKfunct is evaluated in the vector common angles.}

\item{correction}{A character vector specifying the edge correction(s) to be applied among "border", "modified.border", "translate" and "none". The default is "border".}
}
\value{
A list containing:
\itemize{
\item \code{AKhat}: Array containing \code{nds} by \code{ndt} matrices whose elements are \eqn{\hat{K}_{\phi}(r,t)} evaluated in each ang-vector.
\item \code{ds}: If \code{dist} is missing, a vector of distances \code{u} at which \eqn{\hat{K}_{\phi}(u,v)} is computed.
\item \code{dt}: If \code{times} is missing, a vector of distances \code{v} at which \eqn{\hat{K}_{\phi}(u,v)} is computed.
\item \code{lambda}: Value of the estimation of \eqn{\hat{\rho}^{2}}.
\item \code{ang}: If \code{ang} is missing, the function si evalueten on the vector \eqn{(pi/6,pi/4,pi/3,pi/2,2pi/3,3pi/4,5pi/6,pi)}.
}
}
\description{
This function compute an estimation of the anisotropic spatio-temporal inhomogeneous \eqn{K}-function.
}
\examples{
## Not run:
#################

# Realisations of the homogeneous spatio-temporal Poisson processes
stp <- rpp(100)$xyt
# Generated spatio-temporal point pattern
plot(stp)

# Estimation of the anisotropic homogeneous spatio-temporal K-functions
out <- ASTIKhat(stp)

z1 <- out$AKhat[,,1]
z2 <- out$AKhat[,,2]
z3 <- out$AKhat[,,3]
z4 <- out$AKhat[,,4]

# Plot
#par(mfrow=c(2,2))
#out$ds,out$dt,z1,theta=-45,phi=30,zlim=range(z1,na.rm=TRUE),
#ticktype="detailed",xlab="\\n r = #distance",ylab="\\n t = time",
#zlab="",nticks=6,cex.axis=1.3,cex.lab=1.7

## End(Not run)
#################
}
\author{
Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
}
\references{
Gabriel, E., Rowlingson, B., Diggle P J. (2013). \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.

Illian, J. B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Gonzalez, J. A., Rodriguez-Cortes, F. J., Cronie, O., Mateu, J. (2016). Spatio-temporal point process statistics: a review. Spatial Statistics. Accepted.

Ohser, J. and D. Stoyan (1981). On the second-order and orientation analysis of planar stationary point processes. Biometrical Journal 23, 523-533.
}

