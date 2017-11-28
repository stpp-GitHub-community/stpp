\name{rinter}
\alias{rinter}
\title{Generate interaction point patterns}

\description{
Generate one (or several) realisation(s) of the inhibition or contagious process 
in a region \eqn{S\times T}{S x T}.
}

\usage{
 rinter(npoints,s.region,t.region,hs="step",gs="min",thetas=0,
 deltas,ht="step",gt="min",thetat=1,deltat,recent="all",nsim=1,
 discrete.time=FALSE,replace=FALSE,inhibition=TRUE,...)
}

\arguments{
  \item{npoints}{Number of points to simulate. }
  \item{s.region}{Two-column matrix specifying polygonal region containing
  all data locations. 
   If \code{s.region} is missing, the unit square is considered.}
  \item{t.region}{Vector containing the minimum and maximum values of
  the time interval. 
   If \code{t.region} is missing, the interval \eqn{[0,1]}{[0,1]} is considered.} 
  \item{hs, ht}{Function which depends on the distance between points
  and \code{theta}. Can be chosen among "\code{step}" and "\code{gaussian}" 
  or can refer to a user defined function which only depend on d, theta, and delta
  (see details). 
  If \code{inhibition=TRUE}, \code{h} is monotone, increasing, and must tend
  to 1 when the distance tends to infinity. 0 \eqn{\leq}{<=}\code{h(d,theta)}\eqn{\leq}{<=} 1. 
   Otherwise, \code{h} is monotone, decreasing, and must tend
  to 1 when the distance tends to 0.}
  \item{thetas, thetat}{Parameters of \code{hs} and \code{ht} functions.} 
  \item{deltas, deltat}{Spatial and temporal distance of inhibition.}
  \item{gs, gt}{Compute the probability of acceptance of a new point from
  \code{hs} or \code{ht} and \code{recent}. Must be choosen among "\code{min}",
   "\code{max}" and "\code{prod}". }
  \item{recent}{If ``\code{all}'' consider all previous events. If
  is an integer, say \eqn{N}{N}, consider only  the \eqn{N}{N} most recent events.}
  \item{nsim}{Number of simulations to generate. Default is 1. }
  \item{discrete.time}{If TRUE, times belong to \eqn{N}{N},
  otherwise belong to \eqn{{\bf R}^+}{R^+}.} 
  \item{replace}{Logical. If TRUE allows times repeat.}
  \item{inhibition}{Logical. If TRUE, an inhibition process is
  generated. Otherwise, it is a contagious process. }
  \item{...}{Additional parameters if \code{hs} and \code{ht} are defined by the user.}
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
Edith Gabriel <edith.gabriel@univ-avignon.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{plot.stpp}}, \code{\link{animation}} and \code{\link{stan}} for plotting space-time point patterns.
 }
 
\examples{
# simple inhibition process
inh1 = rinter(npoints=200,thetas=0,deltas=0.05,thetat=0,deltat=0.001,
inhibition=TRUE)
plot(inh1$xyt,style="elegant")

\donttest{
# inhibition process using hs and ht defined by the user
hs = function(d,theta,delta,mus=0.1)
{
 res=NULL
 a=(1-theta)/mus
 b=theta-a*delta
 for(i in 1:length(d))
	{	
	if (d[i]<=delta) res=c(res,theta)
	if (d[i]>(delta+mus)) res=c(res,1)
	if (d[i]>delta & d[i]<=(delta+mus)) res=c(res,a*d[i]+b)
	}
 return(res)
}
ht = function(d,theta,delta,mut=0.3)
{
 res=NULL
 a=(1-theta)/mut
 b=theta-a*delta
 for(i in 1:length(d))
	{	
	if (d[i]<=delta) res=c(res,theta)
	if (d[i]>(delta+mut)) res=c(res,1)
	if (d[i]>delta & d[i]<=(delta+mut)) res=c(res,a*d[i]+b)
	}
 return(res)
}
d=seq(0,1,length=100)
plot(d,hs(d,theta=0.2,delta=0.1,mus=0.1),xlab="",ylab="",type="l",
ylim=c(0,1),lwd=2,las=1)
lines(d,ht(d,theta=0.1,delta=0.05,mut=0.3),col=2,lwd=2)
legend("bottomright",col=1:2,lty=1,lwd=2,legend=c(expression(h[s]),
expression(h[t])),bty="n",cex=2)

inh2 = rinter(npoints=100, hs=hs, gs="min", thetas=0.2, deltas=0.1, 
ht=ht, gt="min", thetat=0.1, deltat=0.05, inhibition=TRUE)
animation(inh2$xyt, runtime=15, cex=0.8)

# simple contagious process for given spatial and temporal regions
data(northcumbria)
cont1 = rinter(npoints=100, s.region=northcumbria, t.region=c(1,200), 
thetas=0, deltas=10000, thetat=0, deltat=10, recent=1, inhibition=FALSE)
plot(cont1$xyt,pch=19,s.region=cont1$s.region,mark=TRUE,mark.col=4)
}
}


