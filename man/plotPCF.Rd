\name{plotPCF} \alias{plotPCF}
\title{Plot the estimation of the Space-Time Inhomogeneous Pair Correlation function}

\description{Contour, image or perspective plot of the Space-Time Inhomogeneous Pair correlation function estimate.}

\usage{plotPCF(PCF,n=15,type="contour",legend=TRUE,which=NULL,
main=NULL,...) }

\arguments{
  \item{PCF}{Result of the \code{PCFhat} function.}
    \item{n}{Number of contour levels desired.}
    \item{type}{Specifies the kind of plot: \code{contour} by default, but can also be \code{persp} or \code{image}}
    \item{legend}{Logical indicating whether a legend must be added to the plot.}
    \item{which}{A character specifying the edge correction among the ones used in \code{PCFhat}. If a single edge correction method was used in \code{PCFhat}, it is not necessary to specify \code{which}.}
\item{main}{Plot title.}
\item{...}{Additional arguments to \code{persp} if \code{persp=TRUE}, such as \code{theta} and \code{phi}.}
 }

\seealso{
 \code{\link{contour}}, \code{\link{persp}}, \code{\link{image}} and \code{\link{PCFhat}} for an example.
 }

\author{
Edith Gabriel <edith.gabriel@inrae.fr> }


