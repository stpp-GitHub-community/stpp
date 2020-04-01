\name{rpcp}
\alias{rpcp}
\title{Generate Poisson cluster point patterns}

\description{
Generate one (or several) realisation(s) of the Poisson cluster process in a region 
\eqn{S\times T}{S x T}.
}

\usage{
 rpcp(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, 
 mc=NULL, nsim=1, cluster="uniform", dispersion, infectious=TRUE, 
 edge = "larger.region", larger.region=larger.region, tronc=1,...)
}

\arguments{
  \item{s.region}{Two-column matrix specifying polygonal region containing
  all data locations.  
   If \code{s.region} is missing, the unit square is considered.} 
  \item{t.region}{Vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing, the interval \eqn{[0,1]}{[0,1]} is 
  considered.}
  \item{nparents}{Number of parents. If \code{NULL}, \code{nparents} is from a
  Poisson distribution with intensity \code{lambda}.} 
  \item{npoints}{Number of points to simulate. If \code{NULL} (default), the
  number of points is from a Poisson distribution with mean the double integral 
  of the intensity over \code{s.region} and \code{t.region}.} 
  \item{lambda}{Intensity of the parent process. Can be either a numeric
  value, a function, or a 3d-array (see \code{\link{rpp}}). 
   If \code{NULL}, it is constant and equal to \code{nparents} / volume
  of the domain.} 
  \item{mc}{Average number of children per parent. It is used when
  \code{npoints} is \code{NULL}.}  
  \item{nsim}{Number of simulations to generate. }
  \item{cluster}{Distribution of children: ``uniform'', ``normal'' and
  ``exponential'' are currently implemented. 
  Either a single value if the distribution in space and time is the 
   same, or a vector of length 2, giving first the spatial distribution of
  children  and then the temporal distribution.}
  \item{dispersion}{Scale parameter. It equals twice the standard deviation of 
  location of children relative to their parent for a normal distribution of 
  children; the mean for an exponential distribution and half range for an uniform 
  distribution.} 
  \item{infectious}{If \code{TRUE}, offspring's times are always greater than 
  parent's time.}
  \item{edge}{Specify the edge correction to use "larger.region" or "without".} 
  \item{larger.region}{By default, the larger spatial region is the convex hull 
  of \code{s.region} enlarged by the spatial related value of \code{dispersion} 
  and the larger time interval is \code{t.region} enlarged by the temporal related 
  value of \code{dispersion}.
    One can over-ride default using the 2-vector parameter \code{larger.region}.}
    \item{tronc}{Parameter of the truncated exponential distribution for the distribution of children.}
  \item{...}{Additional parameters of the intensity of the parent process.} 
}

\value{
A list containing:
\item{xyt}{Matrix (or list of matrices if \code{nsim}>1)
containing the points \eqn{(x,y,t)}{(x,y,t)} of the simulated point pattern.
\code{xyt} (or any element of the list if \code{nsim}>1) is an object 
of the class \code{stpp}.}
\item{s.region, t.region}{Parameters passed in argument.}
}

\author{
Edith Gabriel <edith.gabriel@inrae.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{plot.stpp}}, \code{\link{animation}} and \code{\link{stan}} for plotting space-time point patterns.
 }

\examples{
# homogeneous Poisson distribution of parents

data(northcumbria)
pcp1 <- rpcp(nparents=50, npoints=500, s.region=northcumbria, 
t.region=c(1,365), cluster=c("normal","exponential"), 
maxrad=c(5000,5))
\donttest{
animation(pcp1$xyt, s.region=pcp1$s.region, t.region=pcp1$t.region,
runtime=5)
}
# inhomogeneous Poisson distribution of parents

lbda <- function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
pcp2 <- rpcp(nparents=50, npoints=500, cluster="normal", lambda=lbda, 
a=4000/((1-exp(-4))*(1-exp(-2))))
plot(pcp2$xyt, style="elegant")

}
