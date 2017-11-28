\name{animation}
\alias{animation}
\title{Space-time data animation}

\description{
Provide an animation of spatio-temporal point patterns.}

\usage{
animation(xyt, s.region, t.region, runtime=1, incident="red", 
prevalent="pink3", pch=19, cex=0.5, plot.s.region=TRUE, 
scales=TRUE, border.frac=0.05, add=FALSE)
}

\arguments{
  \item{xyt}{Data-matrix containing the \eqn{(x,y,t)}{(x,y,t)}-coordinates. }
  \item{s.region}{Two-column matrix specifying polygonal region containing
  all data-locations \code{xyt[,1:2]}.
   If missing, \code{s.region} is the bounding box of \code{xyt[,1:2]}.}
  \item{t.region}{Interval containing all data-times \code{xyt[,3]}. 
   If missing, \code{t.region} is defined by the range of \code{xyt[,3]}.} 
  \item{runtime}{Approximate running time of animation, in seconds, 
   but it is longer than expected. Can also be \code{NULL}.}  
  \item{incident}{Colour in which incident point \code{xyt[i,1:2]} is
  plotted at time \code{xyt[i,3]}.} 
  \item{prevalent}{Colour to which prevalent point \code{xyt[i,1:2]} fades
  at time \code{xyt[i+1,3]}.} 
  \item{pch}{Plotting symbol used for each point.}
  \item{cex}{Magnification of plotting symbol relative to standard size.}
  \item{plot.s.region}{If \code{TRUE}, plot \code{s.region} as polygon.} 
  \item{scales}{If \code{TRUE}, plot X and Y axes with scales.} 
  \item{border.frac}{Extent of border of plotting region surounding
  \code{s.region}, as fraction of ranges of X and Y.}
  \item{add}{If \code{TRUE}, add the animation to an existing plot.}
 }

\value{
None
}

\author{Peter J Diggle, Edith Gabriel <edith.gabriel@univ-avignon.fr>.}
