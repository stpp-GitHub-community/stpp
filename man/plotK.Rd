\name{plotK} \alias{plotK}
\title{Plot the estimation of the Space-Time Inhomogeneous K-function}

\description{Contour plot or perspective plot or image of the Space-Time
Inhomogeneous K-function estimate.}

\usage{plotK(K,n=15,L=FALSE,type="contour",legend=TRUE,which=NULL,
main=NULL,...) }

\arguments{
  \item{K}{Result of the \code{STIKhat} function.}
    \item{n}{Number of contour levels desired.}
  \item{L}{Logical indicating whether \eqn{K_{ST}(u,v)}{K(u,v)} or \eqn{L(u,v)=K_{ST}(u,v)-\pi u^2 v}{K(u,v)-pi u^2 v} must be plotted.}
    \item{type}{Specifies the kind of plot: \code{contour} by default, but can also be \code{persp} or \code{image}}
    \item{legend}{Logical indicating whether a legend must be added to the plot.}
    \item{which}{A character specifying the edge correction among the ones used in \code{STIKhat}. If a single edge correction method was used in \code{STIKhat}, it is not necessary to specify \code{which}.}
\item{main}{Plot title.}
\item{...}{Additional arguments to \code{persp} if \code{persp=TRUE}, such as \code{theta} and \code{phi}.}
 }

\seealso{
 \code{\link{contour}}, \code{\link{persp}}, \code{\link{image}} and \code{\link{STIKhat}} for an example.
 }

\author{
Edith Gabriel <edith.gabriel@inrae.fr> }


