\name{plotPCF} \alias{plotPCF}
\title{Plot the estimation of the Space-Time Inhomogeneous pair correlation function}

\description{Contour, image or perspective plot of the Space-Time
Inhomogeneous pair correlation function estimate.}

\usage{plotPCF(PCF,n=15,type="contour",legend=TRUE,which=NULL,
main=NULL,...) }

\arguments{
  \item{PCF}{result of the \code{PCFhat} function.}
    \item{n}{number of contour levels desired.}
    \item{type}{specifies the kind of plot: \code{contour} by default, but can also be \code{persp} or \code{image}}
    \item{legend}{logical indicating whether a legend must be added to the plot.}
    \item{which}{A character specifying the edge correction among the ones used in \code{PCFhat}. If a single edge correction method was used in \code{PCFhat}, it is not necessary to specify \code{which}.}
\item{main}{plot title.}
\item{...}{additional arguments to \code{persp} if \code{persp=TRUE}, such as \code{theta} and \code{phi}.}
 }

\seealso{
 \code{\link{contour}}, \code{\link{persp}}, \code{link{image}} and \code{\link{PCFhat}} for an example.
 }

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> }


