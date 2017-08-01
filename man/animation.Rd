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
  \item{xyt}{data-matrix containing the (x,y,t)-coordinates. }
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data-locations \code{xyt[,1:2]}.
   If missing, \code{s.region} is the bounding box of \code{xyt[,1:2]}.}
  \item{t.region}{interval containing all data-times \code{xyt[,3]}. 
   If missing, \code{t.region} is defined by the range of \code{xyt[,3]}.} 
  \item{runtime}{approximate running time of animation, in seconds, 
   but it is longer than expected. Can also be \code{NULL}.}  
  \item{incident}{colour in which incident point \code{xyt[i,1:2]} is
  plotted at time \code{xyt[i,3]}.} 
  \item{prevalent}{colour to which prevalent point \code{xyt[i,1:2]} fades
  at time \code{xyt[i+1,3]}.} 
  \item{pch}{plotting symbol used for each point.}
  \item{cex}{magnification of plotting symbol relative to standard size.}
  \item{plot.s.region}{if \code{TRUE}, plot \code{s.region} as polygon.} 
  \item{scales}{if \code{TRUE}, plot X and Y axes with scales.} 
  \item{border.frac}{extent of border of plotting region surounding
  \code{s.region}, as fraction of ranges of X and Y.}
  \item{add}{if \code{TRUE}, add the animation to an existing plot.}
 }

\value{
None
}

\author{Peter J Diggle, Edith Gabriel <edith.gabriel@univ-avignon.fr>.}
