\name{rpp}
\alias{rpp}
\title{Generate Poisson point patterns}

\description{
Generate one (or several) realisation(s) of the (homogeneous or inhomogeneous) Poisson process in a region \eqn{S\times T}{S x T}.
}

\usage{
rpp(lambda, s.region, t.region, npoints=NULL, nsim=1, replace=TRUE,
    discrete.time=FALSE, nx=100, ny=100, nt=100, lmax=NULL, ...)
}

\arguments{
  \item{lambda}{Spatio-temporal intensity of the Poisson process.
  If \code{lambda} is a single positive number, the function generates 
  realisations of a homogeneous Poisson process, whilst if \code{lambda} 
  is a function of the form \eqn{\lambda(x,y,t,\dots)}{lambda(x,y,t,...)} 
  or a 3D-array it generates realisations of an inhomogeneous Poisson process.}
  \item{s.region}{Two-column matrix specifying polygonal region containing
  all data locations. 
   If \code{s.region} is missing, the unit square is considered.}
  \item{t.region}{Vector containing the minimum and maximum values of
  the time interval. 
   If \code{t.region} is missing, the interval \eqn{[0,1]}{[0,1]} is considered.}
  \item{replace}{Logical allowing times repeat.} 
  \item{npoints}{Number of points to simulate. If \code{NULL}, the
  number of points is from a 
   Poisson distribution with mean the double integral of \code{lambda} over 
   \code{s.region} and \code{t.region}.} 
  \item{discrete.time}{If TRUE, times belong to \eqn{{\bf N}}{N},
  otherwise belong to \eqn{{\bf R}^+}{R^+}.} 
  \item{nsim}{Number of simulations to generate. Default is 1.}
  \item{nx,ny,nt}{Define the size of the 3-D grid on which the intensity
  is evaluated. }
  \item{lmax}{Upper bound for the value of \eqn{\lambda(x,y,t)}{lambda(x,y,t)}, if
  \code{lambda} is a function.} 
  \item{...}{Additional parameters if \code{lambda} is a function.}
}

\value{
A list containing:
\item{xyt}{Matrix (or list of matrices if \code{nsim}>1)
containing the points \eqn{(x,y,t)} of the simulated point pattern.
\code{xyt} (or any element of the list if \code{nsim}>1) is an object 
of the class \code{stpp}.}
\item{Lambda}{\eqn{nx \times ny \times nt}{nx x ny x nt} array of the intensity surface at each time.}
\item{s.region, t.region, lambda}{parameters passed in argument.}
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> and Peter J Diggle.
}

\seealso{
 \code{\link{plot.stpp}}, \code{\link{animation}} and \code{\link{stan}} for plotting space-time point patterns.
 }

\examples{
# Homogeneous Poisson process
# ---------------------------
hpp1 <- rpp(lambda=200,replace=FALSE)

\donttest{stan(hpp1$xyt)}

# fixed number of points, discrete time, with time repeat.
data(northcumbria)
hpp2 <- rpp(npoints=500, s.region=northcumbria, t.region=c(1,1000), 
discrete.time=TRUE)
plot(hpp2$xyt, style="elegant")
\donttest{
polymap(northcumbria)
animation(hpp2$xyt, s.region=hpp2$s.region, t.region=hpp2$t.region, 
runtime=10, add=TRUE)
}

\donttest{
# Inhomogeneous Poisson process
# -----------------------------

# intensity defined by a function
lbda1 = function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
ipp1 = rpp(lambda=lbda1, npoints=400, a=3200/((1-exp(-4))*(1-exp(-2))))
stan(ipp1$xyt)

# intensity defined by a matrix
data(fmd)
data(northcumbria)
h = mse2d(as.points(fmd[,1:2]), northcumbria, nsmse=30, range=3000)
h = h$h[which.min(h$mse)]
Ls = kernel2d(as.points(fmd[,1:2]), northcumbria, h, nx=100, ny=100)
Lt = dim(fmd)[1]*density(fmd[,3], n=200)$y
Lst=array(0,dim=c(100,100,200))
for(k in 1:200) Lst[,,k] <- Ls$z*Lt[k]/dim(fmd)[1]
ipp2 = rpp(lambda=Lst, s.region=northcumbria, t.region=c(1,200), 
discrete.time=TRUE)
           
par(mfrow=c(1,1))
image(Ls$x, Ls$y, Ls$z, col=grey((1000:1)/1000)); polygon(northcumbria)
animation(ipp2$xyt, add=TRUE, cex=0.5, runtime=15)
}
}
