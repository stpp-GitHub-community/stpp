\name{ASTIKhat}
\alias{ASTIKhat}
\title{Anisotropic Space-Time Inhomogeneous \eqn{K}{K}-function}
\usage{
ASTIKhat(xyt, s.region, t.region, lambda, dist, times, ang,
  correction = "border")
}
\arguments{
\item{xyt}{Coordinates and times \eqn{(x,y,t)}{(x,y,t)} of the point pattern.}

\item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}

\item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}

\item{dist}{Vector of distances \eqn{u}{u} at which \eqn{\widehat{K}_{\phi}(r,t)}{K_\phi(r,t)} is computed. If missing, the maximum of \code{dist} is given by \eqn{\min(S_x,S_y)/4}{min(S_x, S_y)/4}, where \eqn{S_x}{S_x} and \eqn{S_y}{S_y} represent the maximum width and height of the bounding box of \code{s.region}.}

\item{times}{Vector of times \eqn{v}{v} at which \eqn{\widehat{K}_{\phi}(r,t)}{K_\phi(r,t)} is computed. If missing, the maximum of \code{times} is given by \eqn{(T_{\max} - T_{\min})/4}{(T_max - T_min)/4}, where \eqn{T_{\min}}{T_min} and \eqn{T_{\max}}{T_max} are the minimum and maximum of the time interval \eqn{T}{T}.}

\item{lambda}{Vector of values of the space-time intensity function evaluated at the points \eqn{(x,y,t)}{(x,y,t)} in \eqn{S\times T}{S x T}. If \code{lambda} is missing, the estimate of the anisotropic space-time \eqn{K}{K}-function is computed as for the homogeneous case, i.e. considering \eqn{n/|S\times T|}{n/|S x T|} as an estimate of the space-time intensity.}

\item{ang}{Angle in radians at which \eqn{\widehat{K}_{\phi}(r,t)}{K_\phi(r,t)} is computed. The argument \code{ang=2*pi} by default.}

\item{correction}{A character vector specifying the edge correction(s) to be applied among "border", "modified.border", "translate" and "none" (see \code{\link{STIKhat}}). The default is "border".}
}
\value{
A list containing:

\item{AKhat}{\code{ndist} x \code{ntimes} matrix containing values of \eqn{\widehat{K}_{\phi}(u,t)}{K_\phi(u,v)}.}
\item{dist, times}{Parameters passed in argument.}
\item{correction}{The name(s) of the edge correction method(s) passed in argument.}
}

\description{
Compute an estimation of the Anisotropic Space-Time inhomogeneous \eqn{K}{K}-function.
}

\author{
Francisco J. Rodriguez-Cortes <frrodriguezc@unal.edu.co> 
}
\references{
Illian, J. B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Gonzalez, J. A., Rodriguez-Cortes, F. J., Cronie, O., Mateu, J. (2016). Spatio-temporal point process statistics: a review. Spatial Statistics. Accepted.

Ohser, J. and D. Stoyan (1981). On the second-order and orientation analysis of planar stationary point processes. Biometrical Journal 23, 523-533.
}