\name{stpp}
\docType{package}
\alias{stpp-package}
\alias{stpp}

\title{Space-Time Point Pattern simulation, visualisation and analysis}

\description{ This package provides models of spatio-temporal point processes in a region \eqn{S\times T}{S x T} and statistical tools for analysing second-order properties of such processes. It also includes static and dynamic (2D and 3D) plots. \code{stpp} is the first dedicated unified computational environment in the area of spatio-temporal point processes.

The \code{stpp} package depends upon some other packages:

\code{splancs}: spatial and space-time point pattern analysis

\code{rgl}: interactive 3D plotting of densities and surfaces

\code{rpanel}: simple interactive controls for R using \code{tcltk} package

\code{KernSmooth}: functions for kernel smoothing for Wand & Jones (1995)

\code{plot3D}: Tools for plotting 3-D and 2-D data
}

\details{
\code{stpp} is a package for simulating, analysing and visualising patterns of points in space and time.

Following is a summary of the main functions and the dataset in the \code{stpp} package.

\emph{To visualise a spatio-temporal point pattern}

\itemize{
\item \code{\link{animation}}: space-time data animation.
\item \code{\link{as.3dpoints}}: create data in spatio-temporal point format.
\item \code{\link{plot.stpp}}: plot spatio-temporal point object. Either a two-panels plot showing spatial locations and cumulative times, or a one-panel plot showing spatial locations with times treated as a quantitative mark attached to each location.
\item \code{\link{stan}}: 3D space-time animation.
}

\emph{To simulate spatio-temporal point patterns}

\itemize{
\item \code{\link{rinfec}}: simulate an infection point process,
\item \code{\link{rinter}}: simulate an interaction (inhibition or contagious) point process,
\item \code{\link{rlgcp}}: simulate a log-Gaussian Cox point process,
\item \code{\link{rpcp}}: simulate a Poisson cluster point process,
\item \code{\link{rpp}}: simulate a Poisson point process,
\item \code{\link{stdcpp}}: simulate a double-cluster point process,
\item \code{\link{sthpcpp}}: simulate a hot-spot point process.
}

\emph{To analyse spatio-temporal point patterns}

\itemize{
\item \code{\link{PCFhat}}: space-time inhomogeneous pair correlation function,
\item \code{\link{STIKhat}}: space-time inhomogeneous K-function.
}

\emph{Dataset}

\code{\link{fmd}}: 2001 food-and-mouth epidemic in north Cumbria (UK).
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr>, Barry Rowlingson, Peter J. Diggle and Francisco J. Rodriguez-Cortes}

\references{Baddeley A., Moller J. and Waagepetersen R. (2000). Non- and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54, 329--350.

Chan, G. and Wood A. (1997). An algorithm for simulating stationary Gaussian random fields. Applied Statistics, Algorithm Section, 46, 171--181.

Chan, G. and Wood A. (1999). Simulation of stationary Gaussian vector fields. Statistics and Computing, 9, 265--268.

Diggle P. , Chedwynd A., Haggkvist R. and Morris S. (1995). Second-order analysis of space-time clustering. Statistical Methods in Medical Research, 4, 124--136.

Gabriel E. (2014). Estimating second-order characteristics of inhomogeneous spatio-temporal point processes: influence ofedge correction methods and intensity estimates. Methodology and computing in Applied Probabillity, 16(1).

Gabriel E., Diggle P. (2009). Second-order analysis of inhomogeneous spatio-temporal point process data. Statistica Neerlandica, 63, 43--51.

Gabriel E., Rowlingson B., Diggle P. (2013). stpp: an R package for plotting, simulating and analyzing 
Spatio-Temporal Point Patterns. Journal of Statistical Software, 
53(2), 1--29.

Gneiting T. (2002). Nonseparable, stationary covariance functions for space-time data. Journal of the American Statistical Association, 97, 590--600.

Gonzalez, J. A., Rodriguez-Cortes, F. J., Cronie, O. and Mateu, J. (2016). Spatio-temporal point process statistics: a review. Spatial Statiscts, 18, 505--544.
}
