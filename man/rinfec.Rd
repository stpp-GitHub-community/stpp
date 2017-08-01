\name{rinfec}
\alias{rinfec}
\title{Generate infection point patterns}

\description{
Generate one (or several) realisation(s) of the infection process 
in a region S x T.
}

\usage{
rinfec(npoints, s.region, t.region, nsim=1, alpha, beta, gamma, 
s.distr="exponential", t.distr="uniform", maxrad, delta, h="step", 
g="min", recent=1, lambda=NULL, lmax=NULL, nx=100, ny=100, nt=1000, 
t0, inhibition=FALSE, ...)
}

\arguments{
  \item{npoints}{number of points to simulate. }
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations. If \code{s.region} is missing, the unit square is considered.} 
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing, the interval [0,1] is considered.} 
  \item{nsim}{number of simulations to generate. Default is 1. }
  \item{alpha}{numerical value for the latent period.} 
  \item{beta}{numerical value for the maximum infection rate.}
  \item{gamma}{numerical value for the infection period.} 
  \item{h}{infection rate function which depends on alpha, beta and delta. 
     Must be choosen among "step" and "gaussian".}
  \item{s.distr}{spatial distribution. Must be choosen among "uniform",
   "gaussian", "exponential" and "poisson". }
  \item{t.distr}{temporal distribution. Must be choosen among "uniform" and
   "exponential".}
  \item{maxrad}{single value or 2-vector of spatial and temporal maximum
  radiation respectively. If single value, the same value is used for space and time.}
   \item{delta}{spatial distance of inhibition/contagion. If missing, the spatial
    radiation is used} 
  \item{g}{Compute the probability of acceptance of a new point from \code{h}
  and \code{recent}. Must be choosen among "min", "max" and "prod". }
  \item{recent}{If ``\code{all}'' consider all previous events. If
  is an integer, say N, consider only the N most recent events.}
  \item{lambda}{function or matrix defining the intensity of a Poisson process if 
  s.distr is Poisson.}
  \item{lmax}{upper bound for the value of lambda.}
  \item{nx, ny}{define the 2-D grid on which the intensity is evaluated if 
  \code{s.distr} is Poisson.} 
  \item{nt}{used to discretize time to compute the infection rate function.}
  \item{t0}{minimum time used to compute the infection rate function. 
  Default is the minimum of \code{t.region}.}
  \item{inhibition}{Logical. If TRUE, an inhibition process is generated. 
   Otherwise, it is a contagious process.}
  \item{...}{additional parameters if \code{lambda} is a function.}
}

\value{
A list containing:
\item{xyt}{matrix (or list of matrices if \code{nsim}>1)
containing the points (x,y,t) of the simulated point pattern.
\code{xyt} (or any element of the list if \code{nsim}>1) is an object 
of the class \code{stpp}.}
\item{s.region, t.region}{parameters passed in argument.}
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{plot.stpp}}, \code{\link{animation}} and \code{\link{stan}} for plotting space-time point patterns.
 }
 
\examples{
\dontrun{
# inhibition; spatial distribution: uniform
inf1 = rinfec(npoints=100, alpha=0.2, beta=0.6, gamma=0.5, 
maxrad=c(0.075,0.5), t.region=c(0,50), s.distr="uniform", 
t.distr="uniform", h="step", p="min", recent="all", 
inhibition=TRUE)
animation(inf1$xyt, cex=0.8, runtime=10)

# contagion; spatial distribution: Poisson with intensity a given matrix
data(fmd)
data(northcumbria)
h = mse2d(as.points(fmd[,1:2]), northcumbria, nsmse=30, range=3000)
h = h$h[which.min(h$mse)]
Ls = kernel2d(as.points(fmd[,1:2]), northcumbria, h, nx=50, ny=50)
inf2 = rinfec(npoints=100, alpha=4, beta=0.6, gamma=20, maxrad=c(12000,20), 
s.region=northcumbria, t.region=c(1,2000), s.distr="poisson", 
t.distr="uniform", h="step", p="min", recent=1, 
lambda=Ls$z, inhibition=FALSE)


image(Ls$x, Ls$y, Ls$z, col=grey((1000:1)/1000)); polygon(northcumbria,lwd=2)
animation(inf2$xyt, add=TRUE, cex=0.7, runtime=15)
}

}
