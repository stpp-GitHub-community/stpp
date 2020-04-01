\name{plot.stpp}
\alias{plot.stpp}
\alias{plot}

\title{Plot for spatio-temporal point objects}

\description{
This function plot either \eqn{xy}{xy}-locations and cumulative distribution of the times, or a space-time 3D scatter, or the time-mark and space-mark of the spatio-temporal point pattern, through arguments \code{style} and \code{type}.

It can also plot \eqn{xy}{xy}-locations with time treated as a quantitative mark attached to each location, as in the previous version of the function, through argument \code{mark} (see stpp version < 2.0.0).
}

\usage{\method{plot}{stpp}(x, s.region=NULL, t.region=NULL, style="generic", type="projection",
mark=NULL , mark.cexmin=0.4, mark.cexmax=1.2, mark.col=1, ...)} 
\arguments{
  \item{x}{Any object of class \code{stpp} in spatio-temporal point format.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If \code{s.region} is missing, the default limits are considered.}
  \item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing,  the default limits are considered.}
  \item{type}{Specify the kind of graphical representation. If \code{type="projection"} (default) the function plot the xy-locations and cumulative distribution of the times. If \code{type="mark"} the function plot the time-mark and space-mark. If \code{type="scatter"} the function plot space-time 3D scatter.}
  \item{style}{Two different classes of graphic styles can be chosen. If \code{style="generic"} (default) the graphics are plot by default function plot in \code{R} and if \code{type="elegant"} the  graphics are plot based on the \code{R} packages ggplot2 and plot3D. }
  \item{mark}{Logical. If \code{NULL} (default), xy-locations and cumulative distribution of the times are plotted. If \code{TRUE}, the time is treated as a quantitative mark attached to each location, and the locations are plotted with the size and/or colour of the plotting symbol determined by the value of the mark.}
  \item{mark.cexmin, mark.cexmax}{Range of the size of the plotting symbol when \code{mark=TRUE}.}
  \item{mark.col}{Colour of the plotting symbol when \code{mark=TRUE}. If \code{mark.col=0}, all locations have the same colour specified by the usual \code{col} argument.
  Otherwise, can be 1 or "black" (default), 2 or "red", 3 or "green", 4 or "blue", in which cases symbols colour is faded, and the darker corresponds to the most recent time.}
  \item{...}{Further arguments to be passed to the functions plot and scatter3D. Typical arguments are \code{pch}, \code{theta} and \code{phi}.}
 }

\value{None}

\seealso{\code{\link{as.3dpoints}} for creating data in spatio-temporal point format.}

\author{Edith Gabriel <edith.gabriel@inrae.fr> and Francisco J. Rodriguez-Cortes.}

\references{
Gabriel E., Rowlingson B., Diggle P. (2013). stpp: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1--29.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. Spatial Statistics. 20, 125-147.}
